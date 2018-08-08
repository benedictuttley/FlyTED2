// <--- HEADER START --->
// Node app Version 1.3
// Status: Logic complete and logging integrated using Winston
// Author: Benedict Uttley
// Date: 04/08/2018
// Current data credit listing: InterMine, FlyMine.
// <--- HEADER END --->

// --- MODULE IMPORTS START --- //
const express = require('express')
const app = express();

const helmet = require('helmet'); // Used to set and ensure HTTP headers correctly for security of application.

const mysql = require('mysql'); // Import mysql module to interact with FlyTed2 database.

const bodyParser = require("body-parser"); // Parse incoming request bodies for example form entry values.

const intermine = require('imjs'); // API to fetch Flymine data.

const handlebars = require('express-handlebars'); // For use in an express environment
const Handlebars = require('handlebars'); // Used as the HTML templating framework to construct the HTML pages server side.

const morgan = require('morgan'); // Used to log http requests.
const winston = require('./config/winston'); // Log various actions and store logs in /logs/app.log

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

app.use(helmet());

app.use(bodyParser.json());

app.use(bodyParser.urlencoded({
  extended: true
}));

app.use(morgan('combined', {
  stream: winston.stream
}));

app.use(express.static('public/img'));
app.use(express.static('public/css'));
app.use(express.static('public/js')); // Gain access to all static files stored in the 'public' directory.
// Serve the 'home' handlebars page when a new user connects to the site, home.hanldebars is the landing page.

app.set('view engine', 'handlebars'); // Create the handlebars engine.

app.listen(3000, () => console.log('App active on port:3000')); // bind application to port 3000;
// --- APP CONFIGURATION END --- //


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

// TODO: Clean below function and if possible modularise.

app.post('/results', (req, res) => { // Listen for incoming probe search requests
  var isEmpty = true;
  FetchExternalData((probe_list) => { // Fetch SQL and FlyMine data
    var asyncTest = function asyncTest(probe, callback) {

      var tissue_query = build_query.createQuery(probe); // Fecth the FlyMine API queries objects

      async.parallel({ // Async module method to enable multiple asynchronous calls all returned in one result object

        // --- STORE THE PROBE START --- //
        probe: callback => { // Also include the probe (CG****) for this result object as a reference move to top.
          callback(null, probe);
        },
        // --- STORE THE PROBE END --- //


        // --- FETCH TISSUE EXPRESSION DATA START --- //
        tissue: callback => {
          flymine.rows(tissue_query).then(rows => {
            if (rows > 0) isEmpty = false;
            callback(null, rows);
          }).catch('error', function(err) {
            //If flymine API error, then log the error:
            winston.error(`${err.status || 500} - ${err.message} - ${req.originalUrl} - ${req.method} - ${req.ip}`);
          });
        },
        // --- FETCH TISSUE EXPRESSION DATA END --- //


        // --- FETCH IMAGE & ANNOTATION DATA START --- //
        flyted: callback => {
          var myQuery = "SELECT * FROM Demo WHERE Probe IN (" + conn.escape(probe) + ")";

          // Check if VARIANT has been set.
          // If yes then Search for entries where the genotypes (a and b)
          // match the users input for the varient.
          if (!(req.body.Variant === undefined || req.body.Variant == "")) {
            myQuery += " AND ( `genotype a` LIKE " + conn.escape('%' + req.body.Variant + '%') +
              " OR `genotype b` LIKE " + conn.escape('%' + req.body.Variant + '%') + ")";
          }

          conn.query(myQuery, (err, my_res) => {
            if (err) {
              //If mysql error, then log the error:
              winston.error(`${err.status || 500} - ${err.message} - ${req.originalUrl} - ${req.method} - ${req.ip}`);
              callback(null, null);
            } else {
              if (my_res.length > 0) isEmpty = false;
              callback(null, my_res);
            }
          });
        },
        // --- FETCH IMAGE & ANNOTATION DATA END --- //

        // --- FETCH THE MICROARRAY ANNOTATIONS AND VALUES START --- //
        affy: function(callback) {
          var myQuery = "SELECT * FROM Probe_Annotations WHERE Gene = (" + conn.escape(probe) + ")";

          conn.query(myQuery, (err, my_res) => {
            if (err) {
              //If mysql error, then log the error:
              winston.error(`${err.status || 500} - ${err.message} - ${req.originalUrl} - ${req.method} - ${req.ip}`);
              callback(null, null);
            } else {
              if (my_res.length > 0) isEmpty = false;
              callback(null, my_res);
            }
          });
        },
        // --- FETCH THE MICROARRAY ANNOTATIONS AND VALUES END --- //

        // --- FETCH THE TRANSCRIPT DATA START --- //
        transcript: function(callback) { // Fetch the probe sequence data from the database.
          var myQuery = "SELECT * FROM Probe_Sequences WHERE Probe IN (" + conn.escape(probe) + ")";

          conn.query(myQuery, (err, my_res) => {
            if (err) {
              //If mysql error, then log the error:
              winston.error(`${err.status || 500} - ${err.message} - ${req.originalUrl} - ${req.method} - ${req.ip}`);
              callback(null, null);
            } else {
              if (my_res.length > 0) isEmpty = false
              callback(null, my_res);
            }
          });
        }
      },
      // --- FETCH THE TRANSCRIPT DATA END --- //

      (err, results) => { // Return result object when all API and SQL functions have executed (all callbacks fired).
        callback(err, results)
      });
    }

    // Async module method to allow async.paralllel method to be placed in for loop so that external data for multiple genes can be required.
    async.map(probe_list, asyncTest, function(err, results) {
      if (err) {
        //If async error, then log the error:
        winston.error(`${err.status || 500} - ${err.message} - ${req.originalUrl} - ${req.method} - ${req.ip}`);
        isEmpty = true; // On error return the 'no results' page, so set isEmpty to true.
      }
      if (isEmpty) res.render('none'); // Return dedicated 'no results' page when no data is found.
      else {
        let data = {
          flyMineData: results, // Stores the API and SQL query results
        }
        res.render('results', data); // Serve Handlebars HTML page.
      }
    });
  }, req.body.Probe, req.body.Variant); // Arguments passed in from user to the async method, these are the keywords of the query.
});

// TODO: Add as an option, perhaps only listing a constant amount each time, [page: queries] ratio?

function FetchExternalData(callback, probes) {
  let probe_list = probes.split(", "); // Create array containing each of the probes entered by the user.
  callback(probe_list);
}


// --- STANDARD ERROR HANDLERS START --- //
// Handle 400 errors.
app.use((req, res) => {
  winston.error(`${404} - ${"FILE NOT FOUND"} - ${req.originalUrl} - ${req.method} - ${req.ip}`);
  res.send("404 FLYTED CANNOT FIND WHAT YOU ARE LOOKING FOR.");
});

// Handle 500 errors.
app.use((err, req, res, next) => {
  winston.error(`${err.status || 500} - ${err.message} - ${req.originalUrl} - ${req.method} - ${req.ip}`);
});

// --- STANDARD ERROR HANDLERS END --- //

// Importing of local modules containing helper functions.
const build_query = require('./queries'); // Contains options for tissue expression query.
const conn = require('./mysql_setup'); // Contains the mysql connection object and authentication details.
