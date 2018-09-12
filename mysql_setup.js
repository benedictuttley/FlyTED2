// --- MODULE DEPENDENCIES START --- //
const morgan = require('morgan'); // Used to log http requests.
const winston = require('./config/winston'); // Log various actions and store logs in /logs/app.log
const mysql = require('mysql'); // Import mysql module to interact with FlyTed2 database.
// --- MODULE DEPENDENCIES END --- //

//--- SQL CONFIGURATION START --- //

let fs = require('fs');
let credentials = JSON.parse(fs.readFileSync('credentials.json','utf8'));

const conn = mysql.createConnection({ // Create mysql connection to FlyTed database.
  host: 'localhost',
  user: credentials.username, // TODO: Change the username to something more meaningful before release.
  password: credentials.password, // TODO: Change the password before release.
  database: 'FlyTed2'
});

// Form and verify connection.
conn.connect(err => {
  if (err) {
    //If mysql connection error, then log the error:
    winston.error(`${err.status || 500} - ${err.message}`);
  } else {
    winston.info("Mysql Connection Established");
  }
});
//--- SQL CONFIGURATION END --- //

module.exports = conn; // Export the connection object.
