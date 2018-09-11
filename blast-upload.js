/* MODULE TO UPLOAD PROBE SEQUENCE DATA AND FIND TARGET SEQUENCE USING PYTHON
   SCRPT THAT UTALISES THE BLAST ALGO: */

const async = require('async');
const mysql = require('mysql');
const winston = require('./config/winston');
const conn = require('./mysql_setup');

module.exports.generate_probe_data = function(filename, callback) {
  var tracker = "";
  var number_successful = 0;
  var unsuccessful_entrys = [];
  var spawn = require('child_process').spawn;
  // Execute the python script:
  var pythonProcesstwo = spawn('python3', ["/home/benedict/test.py", filename]);

  // Listen for python output:
  pythonProcesstwo.stdout.on('data', (data) => {
    console.log("SOME DATA WAS RECEIVED");
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

// Function to upload probe data returned form python script into the Probe_Sequences table:
var blast_mysql = function(row, callback) {
  var values = []
  for (var column in row) {
    if (row.hasOwnProperty(column)) {
      values.push(conn.escape(row[column]));
    }
  }

  let query = "UPDATE Probe_Sequences SET " +
    "`3_Prime_Sequence` = " + conn.escape(row['3-PRIME-SEQUENCE']) +
    ",`5_Prime_Sequence` = " + conn.escape(row['5-PRIME-SEQUENCE']) +
    ",`Target_Sequence_Length` = " + conn.escape(row['TARGET-SEQUENCE-LENGTH']) +
    ",`Target_Sequence` = " + conn.escape(row['TARGET-SEQUENCE']) +
    ",`Transcript_Sequence` = " + conn.escape(row['TRANSCRIPT-SEQUENCE']) +
    ",`Transcript_Sequence_Length` = " + conn.escape(row['TRANSCRIPT-SEQUENCE-LENGTH']) +
    " WHERE Probe = " + conn.escape(row['PROBE']) + "AND Transcript_ID = " + conn.escape(row['TRANSCRIPT']) + ";";

  var myQuery = "INSERT INTO Probe_Sequences(Probe, Transcript_ID, 3_Prime_Sequence," +
    " 5_Prime_Sequence, Target_Sequence_Length, Target_Sequence, Transcript_Sequence," +
    " Transcript_Sequence_Length) VALUES (" + values + ")";

    conn.query(query, (err, res) => {
      // If rows dont exists then create them:
      if (res.affectedRows == 0) {
        conn.query(myQuery, (err, my_res) => {
          if (!err) callback(null, [true, row]);
          else if (err.code == "ER_DUP_ENTRY") {
            winston.error("Duplicate Demo data entered.");
            callback(null, [false, row]);
          } else {
            callback(err, [false, row])
          };
        });
      }

      else{
        if (!err) callback(null, [true, row]);
        else if (err.code == "ER_DUP_ENTRY") {
        winston.error("Duplicate Demo data entered.");
        callback(null, [false, row]);
      } else {
        callback(err, [false, row]);
      }
    }
  });
}
