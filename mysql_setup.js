// --- MODULE DEPENDENCIES START --- //
const morgan = require('morgan'); // Used to log http requests.
const winston = require('./config/winston'); // Log various actions and store logs in /logs/app.log
const mysql = require('mysql'); // Import mysql module to interact with FlyTed2 database.
// --- MODULE DEPENDENCIES END --- //

//--- SQL CONFIGURATION START --- //
const conn = mysql.createConnection({ // Create mysql connection to FlyTed database.
  host: 'localhost',
  user: 'pmauser', // TODO: Change the username to something more meaningful before release.
  password: '[*3G4vlSp7', // TODO: Change the password before release.
  database: 'FlyTed2'
});

// Form and verify connection.
conn.connect(err => {
  if (err) {
    //If mYsql connection error, then log the error:
    winston.error(`${err.status || 500} - ${err.message} - ${req.originalUrl} - ${req.method} - ${req.ip}`);
  } else {
    winston.info("Mysql Connection Established");
  }
});
//--- SQL CONFIGURATION END --- //

module.exports = conn; // Export the connection object.
