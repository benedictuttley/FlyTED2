// Note: Currently the primary key consists of: Slide Name, Genotype A, Genotype B, X-coordinate and
// Y-coordinate.
const express = require('express')
const app = express();
const mysql = require('mysql');
const bodyParser = require("body-parser");

const intermine = require('imjs'); // API to fetch Flymine data.

const handlebars = require('express-handlebars'); // For use in an express environment
const Handlebars = require('handlebars'); // For construncting the HTML files that contain results to user queries.

const async = require('async'); // Used here for multiple asynchronous API calls where all results are stored in one results array.


app.engine('handlebars', handlebars({
  defaultLayout: 'main'
}));
app.set('view engine', 'handlebars');

app.use(bodyParser.json());
app.use(bodyParser.urlencoded({
  extended: true
}));

app.use(express.static('public')); // Provide access to static files stored in the public directory

// Current landing page
app.get('/', (req, res) => {
  res.sendFile('/home/benedict/Desktop/FlyTED2/index.html')
});

// CREDIT: The API queries below were provided through the Flymine site, where UI queries can be translated into their actual query syntax.
// Query [1] --> Find the expression score in all different tissues of the Drosophila for a given Gene:

const flymine = new intermine.Service({
  root: 'http://www.flymine.org/query'
});

app.listen(3000, () => console.log('App active on port:3000'))

// Create mysql connection to existing FlyTed database:
var conn = mysql.createConnection({
  host: 'localhost',
  user: 'pmauser',
  password: '[*3G4vlSp7',
  database: 'FlyTed2'
});

// Test connection:
conn.connect(err => {
  if (err) throw err;
  console.log("Mysql Connection Established");
});

// Query Database for all available rows:
function test(callback) {
  conn.query('SELECT * FROM table_name_for_df ', (err, result) => {
    if (err)
      callback(err, null);
    else
      callback(null, result);
  });
}

// Listen for general query post from HTML form:
app.post('/query', (req, res) => {

  fetchAll((err, result) => {
    if (err) {
      throw err
    } else {
      res.send(result)
    }
  }, req.body.Group);
});



function createExpressionQueries(probeName) {

  let queries = {

    // Expression (tissue) query:
    expression_tissue_query: {
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
        "microArrayResults.tissue.name"
      ],
      "orderBy": [{
        "path": "secondaryIdentifier",
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
    // Expression (stage) query:
    expression_stage_query: {
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

// List for probe query post from HTML form:
app.post('/probe_query', (req, res) => {

  console.log("Probe Query received");

  fetchProbe((err, result, myQuery) => {
    if (err) {
      throw err
    } else {
      JSON.stringify(result);

      // Example async API component A:
      // test: function (callback) {flymine.rows((val) => {callback(null, val);
      var queries = createExpressionQueries(req.body.Probe);

      async.parallel({
        tissue: function(callback) {
          flymine.rows(queries.expression_tissue_query).then(rows => {
            callback(null, rows)
          });
        },
        stage: function(callback) {
          flymine.rows(queries.expression_stage_query).then(rows => {
            callback(null, rows)
          });
        }
      }, function(err, results) {

        let data = {
          query_result: result,
          query: myQuery,
          exp_tissue: results.tissue,
          exp_stage: results.stage
        }
        res.render('hello', data); // Serve Handlebars HTML page.
      });

    }
  }, req.body.Probe);
});

// List for probe query post from HTML form:
app.post('/x_query', (req, res) => {

  console.log("X-coordinate Query received");

  fetchX((err, result) => {
    if (err) {
      throw err
    } else {
      res.send(result)
    }
  }, req.body.x_min, req.body.x_max);
});

// List for probe query post from HTML form:
app.post('/y_query', (req, res) => {

  console.log("Y-coordinate Query received");

  fetchY((err, result) => {
    if (err) {
      throw err
    } else {
      res.send(result)
    }
  }, req.body.yAxis);
});

// List for genotype b query post from HTML form:
app.post('/gena_query', (req, res) => {

  console.log("Genotype A Query received");

  fetchGenA((err, result) => {
    if (err) {
      throw err
    } else {
      res.send(result)
    }
  }, req.body.genA);
});


// List for genotype b query post from HTML form:
app.post('/genb_query', (req, res) => {

  console.log("Genotype B Query received");

  fetchGenB((err, result) => {
    if (err) {
      throw err
    } else {
      res.send(result)
    }
  }, req.body.genB);
});


// TODO: [1] Fetch all entries in table and list them:
function fetchAll(callback, Group) {
  var group_by;
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

// TODO: Fetch all entries and group them by Genotype A:
function group_by_genA() {
  conn.query("SELECT * FROM Demo GROUP BY GenotypeA", (err, result) => {
    if (err) {
      throw err;
    } else {
      callback(null, result);
    }
  });
}

// TODO: Fetch all entries and group them by Genotype B:
function group_by_genB() {
  conn.query("SELECT * FROM Demo GROUP BY GENOTYPE B", (err, result) => {
    if (err) {
      throw err;
    } else {
      callback(null, result);
    }
  });
}


// Select enries corresponding to multiple probes:
function fetchProbe(callback, probeName) {

  let formatted_probes = formatProbes(probeName);
  myQuery = "SELECT * FROM Demo WHERE Probe IN (" + formatted_probes + ")"

  conn.query(myQuery, (err, result) => {
    if (err)
      throw err;
    else {
      callback(null, result, myQuery);
    }
  });
}

function formatProbes(probeName) {
  let formatted_probes = ""
  probeName = probeName.split(", ");
  probeName.forEach(probe => formatted_probes += ("'" + probe + "' "));
  formatted_probes = formatted_probes.trim().replace(/\s+/g, ', ');
  return formatted_probes;
}

// [3] Fetch all entries that have a particular x-coordinate greater than a given number:
function fetchX(callback, xMin, xMax) {

  // Check if min x-value or max x-value input is empty, if true then permit any lower/upper bound:
  xMin.length == 0 ? xMin = 0 : xMin = xMin;
  xMax.length == 0 ? xMax = "(SELECT MAX(Xcoordinate) FROM Demo)" : xMax = xMax;

  conn.query("SELECT * FROM Demo WHERE Xcoordinate >= " + "'" + xMin + "'" + " AND Xcoordinate <= " + xMax, (err, result) => {
    if (err)
      throw err;
    else {
      callback(null, result)
    }
  });
}


// [4] Fetch all entries that have a particular y-coordinate greater than a given number:
function fetchY(callback, yMin, yMax) {

  // Check if min y-value or max y-value input is empty, if true then permit any lower/upper bound:
  yMin.length == 0 ? yMin = 0 : yMin = y_Min;
  yMax.length == 0 ? yMax = "SELECT MAX(Ycoordinate) FROM Demo" : yMax = y_Max;

  conn.query("SELECT * FROM Demo WHERE Ycoordinate >= " + "'" + yMin + "'" + " AND Ycoordinate <= " + yMax, (err, result) => {
    if (err)
      throw err;
    else {
      callback(null, result)
    }
  });
}


// [4] Fetch all entries that have a particular y coordinate greater than a given number:
function fetchGenA(callback, genA) {

  conn.query("SELECT * FROM Demo WHERE GenotypeA = " + "'" + genA + "'", (err, result) => {
    if (err)
      throw err;
    else {
      callback(null, result)
    }
  });
}

// [4] Fetch all entries that have a particular y coordinate greater than a given number:
function fetchGenB(callback, genB) {

  conn.query("SELECT * FROM Demo WHERE GenotypeB = " + "'" + genB + "'", (err, result) => {
    if (err)
      throw err;
    else {
      callback(null, result)
    }
  });
}

// One API to provide limited NCBI data on the gene:
// ncbi.search('gene', 'CG8564', (data) => {console.log(data)});


// Compound query builder:
// Variables to consider: probe, slide, Genotype A and B, Date, Coordinates.

// // TODO: Coordinates Options
// --------------------

// Options for coordinates: The range, from n to p, n to max, as a slider:
// If unchecked then coordinate can be any value

// For FlyMine API: Fetch the following data can be fetched
// Snip location by calculating primer location through performing a BLAST on the 3' and 5'.
// Description
// Expression in different tissues
// Expression at different developmental stages
// Include link to FlyMine
// Include link to FlyBase
