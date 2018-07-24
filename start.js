// <--- HEADER START --->
// Node app Version 1.0
// Status: Incomplete logic but partial core backend functionality implemented
// Author: Benedict Uttley
// Date: 13/07/2018
// Current data credit listing: InterMine, FlyMine, FlyBase, FlyAtlas.
// <--- HEADER END --->

// --- NODE APP SETUP AND MODULE IMPORTS --- //
// Note: Currently the primary key consists of: Slide Name, Genotype A, Genotype B, X-coordinate and Y-coordinate.
const express = require('express')
const app = express();

const mysql = require('mysql'); // Import mysql module to interact with FlyTed2 database.

const bodyParser = require("body-parser"); // Parse incoming request bodies for example form entry values.

const intermine = require('imjs'); // API to fetch Flymine data.

// NOTE: PossiblE FlyMine data to integrate:
// --> Snip location by calculating primer location through performing a BLAST on the 3' and 5'.
// --> Description
// --> Expression in different tissues
// --> Expression at different developmental stages
// --> Include link to FlyMine
// --> Include link to FlyBase

const handlebars = require('express-handlebars'); // For use in an express environment
const Handlebars = require('handlebars'); // Used as the HTML templating framework to construct the HTML pages server side.

const async = require('async'); // Required for multiple asynchronous FlyMine API and mySQL calls where all results are stored in one results array.

//  NOTE: -- CREDIT: The API queries below were provided through the Flymine site,
//  where the respective javascript code for a given user query is generated.

const flymine = new intermine.Service({ // Require the flymine API to integrate data from multiple sources including FlyBase and FlyAtlas.
  root: 'http://www.flymine.org/query' // WARNING: May be moved to HTTPS in the near future!
});

app.engine('handlebars', handlebars({
  defaultLayout: 'main'
}));

app.use(bodyParser.json());

app.use(bodyParser.urlencoded({
  extended: true
}));

app.use(express.static('public')); // Gain access to all static files stored in the 'public' directory.

app.set('view engine', 'handlebars'); // Create the handlebars engine.

// Current landing page
app.get('/', (req, res) => {
  res.sendFile('/home/benedict/Desktop/FlyTED2/index.html');
});

app.listen(3000, () => console.log('App active on port:3000')); // bind application to port 3000;
// --- NODE APP SETUP AND MODULE IMPORTS --- //


// Create mysql connection to FlyTed database.
var conn = mysql.createConnection({
  host: 'localhost',
  user: 'pmauser', // TODO: Change the username to something more meaningful before release.
  password: '[*3G4vlSp7', // TODO: Change the password before release.
  database: 'FlyTed2'
});

// Form connection
conn.connect(err => {
  if (err) throw err;
  console.log("Mysql Connection Established");
});

// Listen for general query post from HTML form:
app.post('/query', (req, res) => {

  fetchAll((err, result) => {
    if (err) {
      throw err;
    } else {
      res.send(result);
    }
  }, req.body.Group);
});

// TODO: Export queries and related helper funnctions such as any formatters to their own module
function createExpressionQueries(probeName) {

  let queries = {

    expression_tissue_query: { // Tissue expresssion query
      "name": "Gene_FlyAtlas",
      "title": "Gene --> FlyAtlas expression data",
      "comment": "06.06.07 updated to work from gene class - Philip",
      "description": "For a given D. melanogaster gene, show expression data from FlyAtlas.",
      "constraintLogic": "A and B and C and D and E and F",
      "from": "Gene",
      "select": [
        "secondaryIdentifier",
        "symbol",
        "microArrayResults.mRNASignal",
        "microArrayResults.mRNASignalSEM",
        "microArrayResults.presentCall",
        "microArrayResults.enrichment",
        "microArrayResults.affyCall",
        "microArrayResults.dataSets.name",
        "microArrayResults.tissue.name",
        "primaryIdentifier"
      ],
      "orderBy": [{
        "path": "microArrayResults.tissue.name",
        "direction": "ASC"
      }],
      "where": [{
          "path": "organism.name",
          "op": "=",
          "value": "Drosophila melanogaster",
          "code": "B",
          "editable": false,
          "switched": "LOCKED",
          "switchable": false
        },
        {
          "path": "microArrayResults",
          "type": "FlyAtlasResult"
        },
        {
          "path": "Gene",
          "op": "LOOKUP",
          "value": probeName,
          "code": "A",
          "editable": true,
          "switched": "LOCKED",
          "switchable": false
        },
        {
          "path": "microArrayResults.affyCall",
          "op": "NONE OF",
          "values": [
            "None"
          ],
          "code": "C",
          "editable": true,
          "switched": "ON",
          "switchable": true
        }
      ]
    },

    expression_stage_query: { // Developmental stage expression query
      "from": "Gene",
      "select": [
        "primaryIdentifier",
        "symbol",
        "rnaSeqResults.stage",
        "rnaSeqResults.expressionScore",
        "rnaSeqResults.expressionLevel"
      ],
      "orderBy": [{
        "path": "rnaSeqResults.stage",
        "direction": "ASC"
      }],
      "where": [{
        "path": "Gene",
        "op": "LOOKUP",
        "value": probeName,
        "extraValue": "D. melanogaster",
        "code": "A",
        "editable": true,
        "switched": "LOCKED",
        "switchable": false
      }]
    }
  }
  return queries;
}

app.post('/probe_query', (req, res) => { // Listen for incoming probe search requests

  console.log("Probe Query received");
  var isEmpty = false;

  FetchExternalData((probe_list, queries) => { // Fetch SQL and FlyMine data
    var asyncTest = function asyncTest(probe, callback) {
      var queries = createExpressionQueries(probe); // Fecth the FlyMine API queries objects
      async.parallel({ // Async module method to enable multiple asynchronous calls all returned in one result object
        tissue: function(callback) { // Fetch tissue expression data for the associated gene.
          flymine.rows(queries.expression_tissue_query).then(rows => {
            callback(null, rows);
          });
        },
        probe: function(callback) { // Also include the probe (CG****) for this result object as a reference
          callback(null, probe);
        },
        flyted: function(callback) { // Fetch images and annotation for the given gene from the mySQL database.
          // TODO: Move SQL logic to dedicated functions and possibly own module.
          var myQuery = "SELECT * FROM Demo WHERE Probe IN (" + ("'" + probe + "'") + ")";

          // Check if VARIANT has been set.
          if (!(req.body.Variant === undefined || req.body.Variant == "")) {
            myQuery += " AND `genotype a` LIKE " + conn.escape('%' + req.body.Variant + '%'); // Check Genotype A column.
            myQuery += " AND `genotype b` LIKE " + conn.escape('%' + req.body.Variant + '%'); // Check Genotype B column.
          }
          conn.query(myQuery, (err, my_res) => {
            if(my_res.length == 0) {
              isEmpty = true;
              console.log("RESULTS ARE VOID");
            }
            if (err) throw err;
            else callback(null, my_res);
          });
        }
      }, (err, results) => {
        callback(null, results)
      });
    }

    // Async module method to allow async.paralllel method to be place in for loop so that external data for multiple genes can be required.
    async.map(probe_list, asyncTest, function(err, results) {
      let data = {
        flyMineData: results, // Stores the API and SQL query results
      }
      console.log(isEmpty);
      if (isEmpty) {
        res.render('none')
      } else {
        res.render('hello', data); // Serve Handlebars HTML page. TODO: Change the name of this handlebars page (e.g queries).
      }
  });
  }, req.body.Probe, req.body.Variant);
});



// TODO: Add this as an option, perhaps only listing a constant amount each time, [page: queries] ratio?
function fetchAll(callback, Group) {
  let group_by = "none";
  Group.length == "None" ? group_by = (" GROUP BY " + Group) : group_by = "";

  conn.query("SELECT * FROM Demo" + group_by, (err, result) => {
    if (err) {
      throw err;
    } else {
      callback(null, result);
    }
  });
}

// TODO: [1] Fetch all entries and group them by probe name:
function group_by_probe() {
  conn.query("SELECT * FROM Demo GROUP BY Probe", (err, result) => {
    if (err) {
      throw err;
    } else {
      callback(null, result);
    }
  });
}

function FetchExternalData(callback, probes) {
  console.log("VVVVVVVVVVVVVVVV");
  let probe_list = probes.split(", "); // Create array containing each of the probes entered by the user onto the form.
  let queries = createExpressionQueries(probe_list); // Fecth the FlyMine API queries objects
  callback(probe_list, queries);
}
