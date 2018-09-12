// <--- HEADER START --->
// Node app Version 1.0
// Status: Code clean up has been performed.
// Author: Benedict Uttley
// Date: 04/08/2018
// Current data credit listing: InterMine, FlyMine.
// <--- HEADER END --->

// --- EXTERNAL MODULE IMPORTS START --- //
const express = require('express')
const app = module.exports = express();

// Admin authentication modules:
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

// Application logger modules:
const morgan = require('morgan'); // Used to log http requests.
const winston = require('./config/winston'); // Custom logging

const async = require('async'); // Required for multiple asynchronous FlyMine API and mySQL calls where all results are stored in one results array.
const formidable = require('formidable');
//  NOTE: -- CREDIT: The tissue expression API queries below were provided through the Flymine site,
//  where the respective javascript code for a given user query is generated.

// File system editor modules:
const fs = require('fs');
const mkdirp = require('mkdirp');
const multer = require('multer');

// Csv writer modules needed as preprocessing to the python blast script:
const createCsvWriter = require('csv-writer').createObjectCsvWriter;

// Require the flymine API to integrate data from multiple sources including FlyBase and FlyAtlas:
const flymine = new intermine.Service({
  root: 'http://www.flymine.org/query' // WARNING: May be moved to HTTPS in the near future!
});
const rc = require('reverse-complement');
// --- EXTERNAL MODULE IMPORTS END --- //


// --- EXTERNAL MODULE CONFIGURATIONS START -- //

var storage = multer.diskStorage({
  destination: function(req, file, cb) {
    cb(null, 'multerTemp')
  },
  filename: function(req, file, cb) {
    cb(null, file.originalname)
  }
})

var testmulter = multer({
  storage: storage
})

// Admin credentials:
const Users = {
  Helen: {
    username: 'Helen',
    password: 'Cardiff_2018'
  }
}

// Hash admin password:
bcrypt.hash(Users.Helen.password, 10, (err, hash) => {
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

app.use(bodyParser.json({
  limit: '150mb'
}));
app.use(bodyParser.urlencoded({ // to support URL-encoded bodies
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

// Permit access to all static files:

app.use(express.static('public/img'));
app.use(express.static('public/css'));
app.use(express.static('public/js'));
app.use(express.static('public/schema'));

app.set('view engine', 'handlebars'); // Create the handlebars engine.


// Define paths for success and faliure upon admin authentication:
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
        if (res) {
          winston.info(`User authentication accepted`);
          return done(null, user);
        } else {
          winston.info(`User authentication failed`);
          return done(null, false, req.flash('error', 'This is a flash message using the express-flash module.'));
        }
      });
    });
  }));

passport.serializeUser((user, done) => { // Create the user for the session.
  var liteUser = createLiteUser(user);
  done(null, liteUser);
});

passport.deserializeUser((user, done) => { // Remove the user from the session.
  done(null, user);
});

function createLiteUser(user) {
  let liteUser = {
    email: "test"
  }
  return liteUser;
}

function findUser(username, done) {
  done(null, Users[username]);
}
// --- EXTERNAL MODULE CONFIGURATIONS END -- //

app.listen(3000, () => winston.info(`FlyTED2 node app is running on port 3000.`)); // bind application to port 3000;

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


// Editor page for admin to make edits to the site.
app.get('/editor', ensureAuthenticated, function(req, res) {
  res.render('editor');
});

// If user is authenticated then dispay the editor page.
app.get('/admin', ensureAuthenticated, (req, res) => {
  res.render('editor');
});


app.post('/fileupload', testmulter.any(), function(req, res) {
  upload.read(req, res, (num_uploaded) => { // Upload submitted image files from admin to correct folder.
    winston.info(`${num_uploaded} image file(s) have been uploaded`);
    if(num_uploaded > 0){
<<<<<<< HEAD
     req.flash('success', `Success: ${num_uploaded} image file(s) have been uploaded`);
    } else{
      req.flash('error', `Error: ${num_uploaded} image file(s) already uploaded.`);
    }

    let flashMessages = res.locals.getMessages();
    return res.status(200).send(flashMessages);
  });
});

app.post('/excel_upload_desc', testmulter.any(), function(req, res) {

  upload.uploadexcel(req, res, (filepath) => {
    annotationuploader.test(conn, filepath, "description", (number_successful, unsuccessful_entrys) => {
      winston.info(`${number_successful} new description entries were successfully inserted.`);

      if (number_successful > 0) req.flash('success', `${number_successful} new image annotation entries were successfully inserted.`);

      unsuccessful_entrys.forEach((entry) => {
        req.flash('error', `The entry for Gene: ${entry['Gene Symbol']} already exists.`);
      });

      let flashMessages = res.locals.getMessages();
      return res.status(200).send(flashMessages);
    });
  });
});

// Listen for bulk image upload and store containing image files in correct directories:
app.post('/excel_upload_image', testmulter.any(), function(req, res) {
  upload.uploadexcel(req, res, (filepath) => {

    annotationuploader.test(conn, filepath, "image_annotations", (number_successful, unsuccessful_entrys) => {
      winston.info(`${number_successful} new image annotation entries were successfully inserted.`);
      if (number_successful > 0) req.flash('success', `${number_successful} new image annotation entries were successfully inserted.`);
      unsuccessful_entrys.forEach((entry) => {
        req.flash('error', `The entry for image with link: [${entry['link_to_file']}] already exists.`);
      });
      let flashMessages = res.locals.getMessages();
      return res.status(200).send(flashMessages);
    });
  });
});

app.post('/excel_probe_data', testmulter.any(), function(req, res) {
  BlastUploader.generate_probe_data(req.files[0].path, (number_successful, unsuccessful_entrys) => {
    winston.info(`${number_successful/2} new probe sequence set(s) were successfully inserted.`);

    if (number_successful > 0) req.flash('success', `${number_successful} new probe sequence entries were successfully inserted.`);
    unsuccessful_entrys.forEach((entry) => {
      req.flash('error', `The entry with file name ${entry['PROBE']} could not be uploaded.`);
});
    let flashMessages = res.locals.getMessages();
    return res.status(200).send(flashMessages);
  });
});

// Helper function for handlebars template creation:
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

// Fetch the microarray data:
app.post('/quick_change_microarray_fetch', (req, res) => {
  let myQuery = "SELECT * FROM Probe_Annotations WHERE Gene = (" + conn.escape(req.body.gene) + ")";

  conn.query(myQuery, (err, my_res) => {
    if (err) winston.error(`${err.status || 500} - ${err.message} - ${req.originalUrl} - ${req.method} - ${req.ip}`);
    return res.status(200).send(my_res);
  });
});

// Fetch the probe data:
app.post('/quick_change_probe_fetch', (req, res) => {

  let query = "SELECT Probe, 3_Prime_Sequence , 5_Prime_Sequence FROM Probe_Sequences WHERE Probe = (" + conn.escape(req.body.gene) + ")";

  conn.query(query, (err, my_res) => {
    if (!err && my_res.length > 0) {
      let sequence = my_res[0]['3_Prime_Sequence'];
      my_res[0]['3_Prime_Sequence'] = rc.reverse_complement(sequence);
      return res.status(200).send(my_res[0]);
    } else {
      if (err) winston.error(`${err.status || 500} - ${err.message} - ${req.originalUrl} - ${req.method} - ${req.ip}`);
      // No resuts sent if error on mysql returned 0 rows:
      return res.status(200).send("");
    }
  });
});

// Update changes to microarray data made by admin:
app.post('/quick_change_microarray_edit', (req, res) => {

  let query = "UPDATE Probe_Annotations SET " +
    "`Probe_Set` = " + conn.escape(req.body['Probe_Set']) +
    ",`Transcript_ID` = " + conn.escape(req.body['Transcript_ID']) +
    ",`Target_Description` = " + conn.escape(req.body['Target_Description']) +
    ",`wt(tin)` = " + conn.escape(req.body['wt(tin)']) +
    ",`wt(mip40)(excised)` = " + conn.escape(req.body['wt(mip40)(excised)']) +
    ",`wt(white)` = " + conn.escape(req.body['wt(mip40)(excised)']) +
    ",`aly5` = " + conn.escape(req.body['aly5']) +
    ",`comr` = " + conn.escape(req.body['comr']) +
    ",`tomb` = " + conn.escape(req.body['tomb']) +
    ",`nxt1` = " + conn.escape(req.body['nxt1']) +
    ",`mip40(sh2)(bsc)` = " + conn.escape(req.body['mip40(sh2)(bsc)']) +
    ",`mip40(ey)` = " + conn.escape(req.body['mip40(ey)']) +
    ",`mip40(ey)(bg4)` = " + conn.escape(req.body['mip40(ey)(bg4)']) +
    ",`wucRNAi` = " + conn.escape(req.body['wucRNAi']) +
    ",`wucRNAi_aly` = " + conn.escape(req.body['wucRNAi_aly']) +
    ",`mip40(sh2)(comr)` = " + conn.escape(req.body['mip40(sh2)(comr)']) +
    ",`mip40(ey)(aly)` = " + conn.escape(req.body['mip40(ey)(aly)']) +
    ",`aly1(27)` = " + conn.escape(req.body['aly1(27)']) +
    ",`aly1(18)` = " + conn.escape(req.body['aly1(18)']) +
    ",`aly1(18)(btr)` = " + conn.escape(req.body['aly1(18)(btr)']) +
    " WHERE Gene = " + conn.escape(req.body['gene']) + ";";

  // Upload new microarray to mysql:
  conn.query(query, (err, my_res) => {

    // If no rows were affected by update query, then insert data as new rows into table:
    if (my_res.affectedRows == 0) {
      var values = [];
      values.push(
        conn.escape(req.body['gene']),
        conn.escape(req.body['Probe_Set']),
        conn.escape(req.body['Transcript_ID']),
        conn.escape(req.body['Target_Description']),
        conn.escape(req.body['wt(tin)']),
        conn.escape(req.body['wt(mip40)(excised)']),
        conn.escape(req.body['wt(white)']),
        conn.escape(req.body['aly5']),
        conn.escape(req.body['comr']),
        conn.escape(req.body['tomb']),
        "0",
        conn.escape(req.body['nxt1']),
        conn.escape(req.body['mip40(sh2)(bsc)']),
        conn.escape(req.body['mip40(ey)']),
        conn.escape(req.body['mip40(ey)(bg4)']),
        conn.escape(req.body['wucRNAi']),
        conn.escape(req.body['wucRNAi_aly']),
        conn.escape(req.body['mip40(sh2)(comr)']),
        conn.escape(req.body['mip40(ey)(aly)']),
        conn.escape(req.body['aly1(27)']),
        conn.escape(req.body['aly1(18)']),
        conn.escape(req.body['aly1(18)(btr)']));

      let aquery = "INSERT INTO Probe_Annotations(Gene,Probe_Set,Transcript_ID," +
        "Target_Description,`wt(tin)`,`wt(mip40)(excised)`,`wt(white)`,aly5,comr," +
        "tomb,nht,nxt1,`mip40(sh2)(bsc)`,`mip40(ey)`,`mip40(ey)(bg4)`,wucRNAi," +
        "wucRNAi_aly,`mip40(sh2)(comr)`,`mip40(ey)(aly)`,`aly1(27)`,`aly1(18)`," +
        "`aly1(18)(btr)`)" + "VALUES (" + values + ")";


      conn.query(aquery, (err, my_res) => {
        if (err) {
          //If mysql error, then log the error:
          req.flash('error', `Quick edit failed: ${err.message}`);
          let flashMessages = res.locals.getMessages();
          winston.error(`${err.status || 500} - ${err.message} - ${req.originalUrl} - ${req.method} - ${req.ip}`);
          return res.status(200).send(flashMessages);
        } else {
          winston.info(`Microarray data has been edited`);
          req.flash('success', `Quick edit was successful.`);
          let flashMessages = res.locals.getMessages();
          return res.status(200).send(flashMessages);
        }
      });

    } else {
      if (err) {
        //If mysql error, then log the error:
        req.flash('error', `Quick edit failed: ${err.message}`);
        let flashMessages = res.locals.getMessages();
        winston.error(`${err.status || 500} - ${err.message} - ${req.originalUrl} - ${req.method} - ${req.ip}`);
        return res.status(200).send(flashMessages);
      } else {
        winston.info(`Microarray data has been edited`);
        req.flash('success', `Quick edit was successful.`);
        let flashMessages = res.locals.getMessages();
        return res.status(200).send(flashMessages);
      }
    }
  });
});

// Update changes to image metadata made by admin:
app.post('/quick_change_img_edit', testmulter.any(), (req, res) => {
  // Upload image files, if present:
  if (req.files[0] != undefined) {
    let image_path = conn.escape(req.body.hidden_original_link);
    upload.copyImage(req.files[0], image_path, () => {
      winston.info(`Image update completed at path: ${image_path}`);
    });
  }

  // Upload new metadata to mysql:

  let query = "UPDATE Image_Data SET " +
    "`file_name` = " + conn.escape(req.body['file_name']) +
    ",`link_to_file` = " + conn.escape(req.body['link_to_file']) +
    ",`slide_name` = " + conn.escape(req.body['slide_name']) +
    ",`date` = " + conn.escape(req.body['date']) +
    ",`user` = " + conn.escape(req.body['user']) +
    ",`probe` = " + conn.escape(req.body['probe']) +
    ",`probe_concentration` = " + conn.escape(req.body['probe_concentration']) +
    ",`genotype_a` = " + conn.escape(req.body['genotype_a']) +
    ",`genotype_b` = " + conn.escape(req.body['genotype_b']) +
    ",`objective` = " + conn.escape(req.body['objective']) +
    ",`optivar` = " + conn.escape(req.body['optivar']) +
    ",`cmount` = " + conn.escape(req.body['cmount']) +
    ",`stages_shown_in_picture` = " + conn.escape(req.body['stages_shown_in_picture']) +
    ",`description_of_staining_pattern` = " + conn.escape(req.body['description_of_staining_pattern']) +
    ",`comments` = " + conn.escape(req.body['comments']) +
    ",`xcoordinate` = " + conn.escape(req.body['xcoordinate']) +
    ",`ycoordinate` = " + conn.escape(req.body['ycoordinate']) +
    ",`Apical_spermatogonia` = " + conn.escape(req.body['Apical_spermatogonia']) +
    ",`Germ_line` = " + conn.escape(req.body['Germ_line']) +
    ",`Somatic` = " + conn.escape(req.body['Somatic']) +
    ",`TARGET_SEQUENCE` = " + conn.escape(req.body['TARGET_SEQUENCE']) +
    " WHERE probe = " + conn.escape(req.body['hidden_original_probe']) + " AND link_to_file = " +
    conn.escape(req.body['hidden_original_link']) + ";";

conn.query(query, (err, my_res) => {

    if (err) {
      //If mysql error, then log the error:
      req.flash('error', `Quick edit failed: ${err.message}`);
      let flashMessages = res.locals.getMessages();
      winston.error(`${err.status || 500} - ${err.message} - ${req.originalUrl} - ${req.method} - ${req.ip}`);
      return res.status(200).send(flashMessages);
    } else {
      req.flash('success', `Quick edit was successful.`);
      let flashMessages = res.locals.getMessages();
      return res.status(200).send(flashMessages);
    }
  });
});

/* When admin uploads probe sequences, add them to the database
after formating and analyis with python script: */

app.post('/quick_change_probe_edit', (req, res) => {

  const csvWriter = createCsvWriter({
    path: 'this_test.csv',
    header: [{
        id: 'probe',
        title: 'PROBE'
      },
      {
        id: 'sequence',
        title: 'SEQUENCE'
      }
    ],
    append: false
  });
  var records = [{
      probe: "3_" + req.body.Probe,
      sequence: req.body['3_Prime_Sequence']
    },
    {
      probe: "5_" + req.body.Probe,
      sequence: req.body['5_Prime_Sequence']
    }
  ];
  csvWriter.writeRecords(records).then(() => {
    BlastUploader.generate_probe_data("this_test.csv", (number_successful, unsuccessful_entrys) => {
      winston.info(`${number_successful} new probe sequence entries were successfully inserted.`);

      if (number_successful > 0) req.flash('success', `${number_successful} probe sequence entries were successfully updated/created.`);
      else req.flash('error', `Probes could not be uploaded.`);

      let flashMessages = res.locals.getMessages();
      return res.status(200).send(flashMessages);
    });
  });
});

// Fetch image metadata to be edited:
app.post('/quick_change_image_fetch', (req, res) => {
  let query = "SELECT * FROM Image_Data WHERE probe = (" + conn.escape(req.body.gene) + ")";
  conn.query(query, (err, my_res) => {
    if (err) winston.error(`${err.status || 500} - ${err.message} - ${req.originalUrl} - ${req.method} - ${req.ip}`);
    return res.status(200).send(my_res);
  });
});

// Add new image form entry to mysql.
app.post('/QuickImageAddition', testmulter.any(), (req, res) => {



  // Add Metadata:
  let row_to_upload = [];
  let row = req.body;
  row_to_upload.push(conn.escape(row['date']), conn.escape(row['Link_To_File']), conn.escape(row['Slide_Name']), conn.escape(row['date']),
    conn.escape(row['user']), conn.escape(row['probe']), conn.escape(row['probe_concentration']), conn.escape(row['genotype_a']),
    conn.escape(row['genotype_b']), conn.escape(row['objective']), conn.escape(row['optivar']), conn.escape(row['cmount']),
    conn.escape(row['stages_shown_in_picture']), conn.escape(row['description_of_staining_pattern']), conn.escape(row['comments']),
    conn.escape(row['xcoordinate']), conn.escape(row['ycoordinate']));

  // Add image to correct path on the file system:
  annotationuploader.quickAdditionImageUpload(row_to_upload, (err) => {
    upload.read(req, res, (num_uploaded) => {
      if (err) {
        winston.error(`${err.status || 500} - ${err.message} - ${req.originalUrl} - ${req.method} - ${req.ip}`);
        req.flash('error', `The addition was bad`)
      } else req.flash('success', `The addition was good`);
      let flashMessages = res.locals.getMessages();
      return res.status(200).send(flashMessages);
    });
  });
});

// Listen for when admin deletes a image/metadata pair:

app.post('/removeImageAndData', (req, res) => {

  //Mysql to delete image and metadata:
  var img_name = conn.escape(req.body.image_name);
  let query = ("DELETE FROM Image_Data WHERE link_to_file = " + img_name);

  conn.query(query, (err, my_res) => {
    if (err) winston.error(`${err.status || 500} - ${err.message} - ${req.originalUrl} - ${req.method} - ${req.ip}`);
    else {
      req.flash('success', `Metadata for image: ${img_name} has been successfully deleted`);
      winston.info(`Metadata for image with file name ${img_name} has been successfully deleted.`);
     }

    // Need to now remove the image:
    let file_to_remove = ("/var/www/html/FlyTED2/public/img" + req.body.image_name);
    fs.unlink(file_to_remove, (err) => {
      if (err) {
        winston.error(`${err.status || 500} - ${err.message} - ${req.originalUrl} - ${req.method} - ${req.ip}`);
        req.flash('error', `The image ${'img_name'} could not be removed.`);
      } else req.flash('success', `The image ${'img_name'} has been sucessfully removed.`);

      if (number_successful > 0) req.flash('success', `${number_successful} probe sequence entries were successfully updated/created.`);
      else req.flash('error', `Probes could not be uploaded.`);

      let flashMessages = res.locals.getMessages();
      return res.status(200).send(flashMessages);
    });
  });
});

// Fetch image metadata to be edited:
app.post('/quick_change_image_fetch', (req, res) => {
  let query = "SELECT * FROM Demo WHERE probe = (" + conn.escape(req.body.gene) + ")";
  conn.query(query, (err, my_res) => {
    if (err) winston.error(`${err.status || 500} - ${err.message} - ${req.originalUrl} - ${req.method} - ${req.ip}`);
    return res.status(200).send(my_res);
  });
});

// Add new image form entry to mysql.
app.post('/QuickImageAddition', testmulter.any(), (req, res) => {



  // Add Metadata:
  let row_to_upload = [];
  let row = req.body;
  row_to_upload.push(conn.escape(row['date']), conn.escape(row['Link_To_File']), conn.escape(row['Slide_Name']), conn.escape(row['date']),
    conn.escape(row['user']), conn.escape(row['probe']), conn.escape(row['probe_concentration']), conn.escape(row['genotype_a']),
    conn.escape(row['genotype_b']), conn.escape(row['objective']), conn.escape(row['optivar']), conn.escape(row['cmount']),
    conn.escape(row['stages_shown_in_picture']), conn.escape(row['description_of_staining_pattern']), conn.escape(row['comments']),
    conn.escape(row['xcoordinate']), conn.escape(row['ycoordinate']));

  // Add image to correct path on the file system:
  annotationuploader.quickAdditionImageUpload(row_to_upload, (err) => {
    upload.read(req, res, (num_uploaded) => {
      if (err) {
        winston.error(`${err.status || 500} - ${err.message} - ${req.originalUrl} - ${req.method} - ${req.ip}`);
        req.flash('error', `The addition was bad`)
      } else req.flash('success', `The addition was good`);
      let flashMessages = res.locals.getMessages();
      return res.status(200).send(flashMessages);
    });
  });
});


app.post('/removeProbeData', (req, res) => {
  var gene_name = conn.escape(req.body.gene);
  let query = ("DELETE FROM Probe_Sequences WHERE Probe = " + gene_name);

conn.query(query, (err, my_res) => {
    if (err) {
      winston.error(`${err.status || 500} - ${err.message} - ${req.originalUrl} - ${req.method} - ${req.ip}`);
      req.flash('error', `The deletion was bad`);
      let flashMessages = res.locals.getMessages();
      return res.status(200).send(flashMessages);
    } else {
      winston.info(`Probe sequence data for gene ${gene_name} has been removed.`);
      req.flash('success', `The deletion was good`);
      let flashMessages = res.locals.getMessages();
      return res.status(200).send(flashMessages);
    }
  });
});

app.post('/removeMicroarrayData', (req, res) => {
  var gene_name = conn.escape(req.body.gene);
  let query = ("DELETE FROM Probe_Annotations WHERE gene = " + gene_name);

  conn.query(query, (err, my_res) => {
    if (err) {
      winston.error(`${err.status || 500} - ${err.message} - ${req.originalUrl} - ${req.method} - ${req.ip}`);
      req.flash('error', `Microarray data for gene ${gene_name} was not successful, Error: ${err.message}`);
    } else {
      winston.info(`Microarray data for gene ${gene_name} has been removed.`);
      req.flash('success', `Microarray data for gene ${gene_name} was deleted successfully`);
    }

    let flashMessages = res.locals.getMessages();
    return res.status(200).send(flashMessages);
  });
});

// --- STANDARD ERROR HANDLERS START ---
// Handle 400 errors.
app.use((req, res) => {
  winston.error(`${404} - ${"FILE NOT FOUND"} - ${req.originalUrl} - ${req.method} - ${req.ip}`);
  res.send("404 FLYTED CANNOT FIND WHAT YOU ARE LOOKING FOR.");
});

// Handle 500 errors.
app.use((err, req, res, next) => {
  winston.error(`${err.status || 500} - ${err.message} - ${req.originalUrl} - ${req.method} - ${req.ip}`);
});
// --- STANDARD ERROR HANDLERS END ---

// Importing of local modules containing helper functions.
const build_query = require('./queries'); // Contains options for tissue expression query.
const conn = require('./mysql_setup'); // Contains the mysql connection object and authentication details.
const upload = require('./upload'); // Local module to save file uploads (currently expression images only);
const annotationuploader = require('./annotation-uploader'); // Upload new annotation data
const BlastUploader = require('./blast-upload'); // module that interacts with py script to perform BLAST analysis on subitted probe sequences

// Check if user is authenticated before serving the editor page:
function ensureAuthenticated(req, res, next) {
  if (req.isAuthenticated()) {
    return next(); // Successful authentication so executed callback
  } else {
    // Access denied, user is redirected back to login:
    req.flash('error', 'You must be logged in to view this resource.');
    res.render('admin', {
      errors: res.locals.getMessages()
    });
  }
}
