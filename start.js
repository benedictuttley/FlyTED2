// <--- HEADER START --->
// Node app Version 1.0
// Status: Code clean up has been performed.
// Author: Benedict Uttley
// Date: 04/08/2018
// Current data credit listing: InterMine, FlyMine.
// <--- HEADER END --->

// --- MODULE IMPORTS START --- //
const express = require('express')
const app = module.exports = express();
const passport = require('passport');
const LocalStrategy = require('passport-local').Strategy;
const expressSession = require('express-session');
const flash = require('express-flash-messages');
const bcrypt = require('bcrypt');
const minifyHTML = require('express-minify-html')
const minify = require('express-minify'); // Minify CSS and JS files.
const helmet = require('helmet'); // Used to set and ensure HTTP headers correctly for security of application.

const mysql = require('mysql'); // Import mysql module to interact with FlyTed2 database.

const bodyParser = require("body-parser"); // Parse incoming request bodies for example form entry values.

const intermine = require('imjs'); // API to fetch Flymine data.

const handlebars = require('express-handlebars'); // For use in an express environment
const Handlebars = require('handlebars'); // Used as the HTML templating framework to construct the HTML pages server side.

const morgan = require('morgan'); // Used to log http requests.
const winston = require('./config/winston'); // Log various actions and store logs in /logs/app.log

const async = require('async'); // Required for multiple asynchronous FlyMine API and mySQL calls where all results are stored in one results array.
const formidable = require('formidable');
//  NOTE: -- CREDIT: The tissue expression API queries below were provided through the Flymine site,
//  where the respective javascript code for a given user query is generated.
const fs = require('fs');
const mkdirp = require('mkdirp');
const multer = require('multer');

var storage = multer.diskStorage({
  destination: function (req, file, cb) {
    cb(null, 'multerTemp')
  },
  filename: function (req, file, cb) {
    cb(null, file.originalname)
  }
})


var testmulter = multer({ storage: storage })
const flymine = new intermine.Service({ // Require the flymine API to integrate data from multiple sources including FlyBase and FlyAtlas.
  root: 'http://www.flymine.org/query' // WARNING: May be moved to HTTPS in the near future!
});
// --- MODULE IMPORTS END --- //

// IN MEMORY USER
const Users = {
  Helen: {
    username: 'Helen',
    password: 'Cardiff_2018'
  }
}

bcrypt.hash(Users.Helen.password, 10, (err, hash) => {
  // Store the hashed password instead.
  Users.Helen.password = hash;
});

// --- APP CONFIGURATION START --- //
app.engine('handlebars', handlebars({
  defaultLayout: 'main'
}));

app.use(helmet());

app.use(minifyHTML({
  override: true,
  exception_url: false,
  htmlMinifier: {
    removeComments: true,
    collapseWhitespace: true,
    collapseBooleanAttributes: true,
    removeAttributeQuotes: true,
    removeEmptyAttributes: true,
    minifyJS: true
  }
}));

app.use(bodyParser.json({limit: '150mb'}));
app.use(bodyParser.urlencoded({     // to support URL-encoded bodies
limit: '150mb',
extended: true
}));
app.use(morgan('combined', {
  stream: winston.stream
}));
app.use(expressSession({
  secret: 'mySecretKey',
  resave: true,
  saveUninitialized: true
}));
app.use(passport.initialize());
app.use(passport.session());
app.use(flash());



app.use(express.static('public/img'));
app.use(express.static('public/css'));
app.use(express.static('public/js')); // Gain access to all static files stored in the 'public' directory.
// Serve the 'home' handlebars page when a new user connects to the site, home.hanldebars is the landing page.

app.set('view engine', 'handlebars'); // Create the handlebars engine.


// TEST PASSPORT ROUTE
app.post('/login', passport.authenticate('login', {
  successRedirect: '/editor',
  failureRedirect: '/admin',
}));


// When user attempts login, authenticate with Passport.js
passport.use("login", new LocalStrategy({
    passReqToCallback: true // Allow Passport.js to use the request object
  },
  (req, username, password, done) => {
    findUser(username, (err, user) => {
      if (err)
        return done(err);
      if (!user) {
        return done(null, false)
      }

      bcrypt.compare(password, user.password, (err, res) => {
        if(res){
          console.log("THE PASSWORDS MATCH!");
          return done(null, user); // On correct credentials
          }
        else{
          console.log("THE PASSWORDS DONT MATCH");
          return done(null, false, req.flash('error', 'This is a flash message using the express-flash module.'));
        }
        });

      });
    }));

passport.serializeUser((user, done) => {  // Create the user for the session.
  var liteUser = createLiteUser(user);
  done(null, liteUser);
});

passport.deserializeUser((user, done) => { // Remove the user from the session.
  done(null, user);
});

function createLiteUser(user){
  let liteUser = {
    email: "test"
  }
  return liteUser;
}

function findUser(username, done) {
  done(null, Users[username]);
}


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

// Serve the 'cite' handlebars page upon request:
app.get('/cite', (req, res) => {
  res.render('cite');
})

// Serve the 'credits' handlebars page upon request:
app.get('/credits', (req, res) => {
  res.render('credits');
})


function test(file) {
  console.log(file.name);
  req.flash('error',file.name);
}

// app.post('/test', function(req, res) {
//   console.log("RECIEVED!");
//   var form = new formidable.IncomingForm();
//   form.multiples = true; // per their documents
//   form.maxFileSize = 200 * 1024 * 1024;
//   form.parse(req, function(err, fields, files) {
//     console.log(files);
//     console.log(err);
//     res.render('home');
//     res.status(201).end()
// });
//
// form.on('file', function(field, file) {
//        console.log("---" + file.name);
//    })
//
//    form.on('error', function(err) {
//      console.log("ihoihoiho");
// });
//
// res.render('home');
// });





// Editor page for admin to make edits to the site.
app.get('/editor', ensureAuthenticated, function(req, res) {
  res.render('editor');
});

// If user is authenticated then dispay the editor page.
app.get('/admin',ensureAuthenticated,(req, res) => {
  res.render('editor');
});


app.post('/fileupload', testmulter.any(), function(req, res) {
  upload.read(req, res, (num_uploaded) => { // Upload submitted image files from admin to correct folder.
    winston.info(`${num_uploaded} image file(s) have been uploaded`);
    let flashMessages = res.locals.getMessages();
    return res.status(200).send(flashMessages);
  });
});


Handlebars.registerHelper('json', function(context) {
  return JSON.stringify(context);
});

app.post('/results', (req, res) => { // Listen for incoming probe search requests
  var isEmpty = true;
  FetchExternalData((probe_list) => { // Fetch SQL and FlyMine data
    var asyncTest = function asyncTest(probe, callback) {

      var tissue_query = build_query.createQuery(probe); // Fecth the FlyMine API queries objects

      async.parallel({ // Async module method to enable multiple asynchronous calls all returned in one result object

          // --- STORE THE PROBE START --- //
          probe: callback => { // Also include the probe (CG****) for this result object as a reference move to top.
            probe = probe.toUpperCase()
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
            var myQuery = "SELECT * FROM Image_Data WHERE Probe IN (" + conn.escape(probe) + ")";

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
          affy: callback => {
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
          transcript: callback => { // Fetch the probe sequence data from the database.
            var myQuery = "SELECT * FROM Probe_Sequences WHERE Probe IN (" + conn.escape(probe) + ")";

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
      if (isEmpty) {
        res.render('none'); // Return dedicated 'no results' page when no data is found.
      } else {
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
const upload = require('./upload'); // Local module to save file uploads (currently expression images only);

function ensureAuthenticated(req, res, next) {
  if (req.isAuthenticated()) {
    return next();
  } else{
    // denied. redirect to login
    req.flash('error', 'You must be logged in to view this resource.');
    var flashMessages = res.locals.getMessages(); // Fetch errors sored in flash memory.
    res.render('admin', {
      errors: flashMessages
    });
  }
}
