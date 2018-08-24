const async = require('async');
const mysql = require('mysql');
const winston = require('./config/winston');
const conn = require('./mysql_setup');

module.exports.generate_probe_data = function(filename, callback) {
  var tracker = "";
  number_successful = 0;
  unsuccessful_entrys = [];

  var spawn = require('child_process').spawn;
  var pythonProcesstwo = spawn('python3', ["/home/benedict/test.py", filename]);

  pythonProcesstwo.stdout.on('data', (data) => {
    tracker += String.fromCharCode.apply(null, data);
  });

  pythonProcesstwo.stdout.on('end', () => {
    number_successful = 0;
    unsuccessful_entrys = [];
    tracker = JSON.parse(tracker);
    async.map(tracker, blast_mysql, function(err, results) {
      results.forEach(result => {
        if (result[0]) number_successful = number_successful + 1;
        else unsuccessful_entrys.push(result[1]);
      });

      callback(number_successful, unsuccessful_entrys);
    });
  });
}

// Function to upload probe data returned form python script into the Probe_Sequences mysql database.
var blast_mysql = function(row, callback) {
  var values = []
  for (var column in row) {
    if (row.hasOwnProperty(column)) {
      values.push("'" + row[column] + "'");
    }
  }

  var myQuery = "INSERT INTO Probe_Sequences(Probe, Transcript_ID, 3_Prime_Sequence," +
    " 5_Prime_Sequence, Target_Sequence_Length, Target_Sequence, Transcript_Sequence," +
    " Transcript_Sequence_Length) VALUES (" + values + ")";

  conn.query(myQuery, (err, res) => {
    if (!err) callback(null, [true, row]);
    else if (err.code == "ER_DUP_ENTRY") {
      winston.error("Duplicate Demo data entered.");
      callback(null, [false, row]);
    } else {
      callback(err, [false, row])
    };
  });
}
