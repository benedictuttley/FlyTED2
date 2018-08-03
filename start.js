// <--- HEADER START --->
// Node app Version 1.2
// Status: Logic complete
// Author: Benedict Uttley
// Date: 02/08/2018
// Current data credit listing: InterMine, FlyMine.
// <--- HEADER END --->

// --- MODULE IMPORTS START --- //
const express = require('express')
const app = express();

const mysql = require('mysql'); // Import mysql module to interact with FlyTed2 database.

const bodyParser = require("body-parser"); // Parse incoming request bodies for example form entry values.

const intermine = require('imjs'); // API to fetch Flymine data.

const handlebars = require('express-handlebars'); // For use in an express environment
const Handlebars = require('handlebars'); // Used as the HTML templating framework to construct the HTML pages server side.

const async = require('async'); // Required for multiple asynchronous FlyMine API and mySQL calls where all results are stored in one results array.

//  NOTE: -- CREDIT: The tissue expression API queries below were provided through the Flymine site,
//  where the respective javascript code for a given user query is generated.

const flymine = new intermine.Service({ // Require the flymine API to integrate data from multiple sources including FlyBase and FlyAtlas.
  root: 'http://www.flymine.org/query' // WARNING: May be moved to HTTPS in the near future!
});
// --- MODULE IMPORTS END --- //

// --- APP CONFIGURATION START --- //
app.engine('handlebars', handlebars({
  defaultLayout: 'main'
}));

app.use(bodyParser.json());

app.use(bodyParser.urlencoded({
  extended: true
}));

app.use(express.static('public')); // Gain access to all static files stored in the 'public' directory.
// Serve the 'home' handlebars page when a new user connects to the site, home.hanldebars is the landing page.

app.set('view engine', 'handlebars'); // Create the handlebars engine.

app.listen(3000, () => console.log('App active on port:3000')); // bind application to port 3000;
// --- APP CONFIGURATION END --- //

// --- SQL CONFIGURATION START --- //
const conn = mysql.createConnection({ // Create mysql connection to FlyTed database.
  host: 'localhost',
  user: 'pmauser', // TODO: Change the username to something more meaningful before release.
  password: '[*3G4vlSp7', // TODO: Change the password before release.
  database: 'FlyTed2'
});

// Form and verify connection.
conn.connect(err => {
  if (err) throw err;
  console.log("Mysql Connection Established");
});
// --- SQL CONFIGURATION END --- //



// MAIN SERVER LOGIC BELOW //

// --- PAGE REQUESTS ---

// Serve the 'home' handlebars page when a new user connects to the site, home.hanldebars is the landing page.
app.get('/', (req, res) => {
  res.render('home');
});

// Serve the 'home' handlebars page upon request:
app.get('/home', (req, res) => {
  res.render('home');
})

// Serve the 'credits' handlebars page upon request:
app.get('/credits', (req, res) => {
  res.render('credits');
})

// Serve the 'cite' handlebars page upon request:
app.get('/cite', (req, res) => {
  res.render('cite');
})

Handlebars.registerHelper('json', function(context) {
  return JSON.stringify(context);
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

app.post('/results', (req, res) => { // Listen for incoming probe search requests
  var isEmpty = false;

  FetchExternalData((probe_list, queries) => { // Fetch SQL and FlyMine data
    var asyncTest = function asyncTest(probe, callback) {
      var queries = createExpressionQueries(probe); // Fecth the FlyMine API queries objects
      async.parallel({ // Async module method to enable multiple asynchronous calls all returned in one result object
        tissue: function(callback) { // Fetch tissue expression data for the associated gene.
          flymine.rows(queries.expression_tissue_query).then(rows => {
            callback(null, rows);
          }).catch('error', function(e) {
            console.log(e); // Check for error with API and log the error.
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
            // Search for entries where the genotypes (a and b) match the users input for the varient.
            myQuery += " AND ( `genotype a` LIKE " + conn.escape('%' + req.body.Variant + '%') +
              " OR `genotype b` LIKE " + conn.escape('%' + req.body.Variant + '%') + ")";
          }

          conn.query(myQuery, (err, my_res) => {
            if (my_res.length == 0) {
              isEmpty = true;
              console.log("RESULTS ARE VOID!");
            }
            if (err) throw err;
            else callback(null, my_res);
          });
        },
        affy: function(callback) { // Fetch the affy data and annotations from the database.
          var myQuery = "SELECT * FROM Probe_Annotations WHERE Gene = (" + ("'" + probe + "'") + ")";

          conn.query(myQuery, (err, my_res) => {
            console.log(my_res);
            if (err) throw err;
            else callback(null, my_res);
          });
        },

        transcript: function(callback) { // Fetch the probe sequence data from the database.
          var myQuery = "SELECT * FROM Probe_Sequences WHERE Probe IN (" + ("'" + probe + "'") + ")";

          conn.query(myQuery, (err, my_res) => {
            if (err) throw err;
            else callback(null, my_res);
          });
        }

      }, (err, results) => { // Return the result of each indivisulal query in one object called results, once all api/db calls have completed.
        callback(null, results)
      });
    }

    // Async module method to allow async.paralllel method to be placed in for loop so that external data for multiple genes can be required.
    async.map(probe_list, asyncTest, function(err, results) {
      let data = {
        flyMineData: results, // Stores the API and SQL query results
      }
      if (isEmpty) {
        res.render('none'); // Return dedicated 'no results' page when no data is found (image-wise) NOTE: Could be adapted for other data sets being void instead.
      } else {
        res.render('results', data); // Serve Handlebars HTML page.
      }
    });
  }, req.body.Probe, req.body.Variant); // Arguments passed in from user to the async method, these are the keywords of the query.
});

// TODO: Add as an option, perhaps only listing a constant amount each time, [page: queries] ratio?

function FetchExternalData(callback, probes) {
  let probe_list = probes.split(", "); // Create array containing each of the probes entered by the user onto the form.
  let queries = createExpressionQueries(probe_list); // Fecth the FlyMine API queries objects
  callback(probe_list, queries);
}

// Handle 404 errors.
app.use(function(req, res) { //TODO: REPLACE WITH CUSTOM ERROR PAGE.
  res.send('404: Page not Found', 404);
});

// Handle 500 errors.
app.use(function(error, req, res, next) { //TODO: REPLACE WITH CUSTOM ERROR PAGE.
  res.send('500: Internal Server Error', 500);
});
