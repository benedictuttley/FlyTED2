//Module to read in annottaions associated with new image uploads.
const xlstojson = require('xls-to-json');
const mysql = require('mysql');
const conn = require('./mysql_setup');
const winston = require('./config/winston'); // Log various actions and store logs in /logs/app.log
const async = require('async');
module.exports.test = function(conn, filename, callback){

xlstojson({
      input: filename, // input xls
      output: "output.json", // output json
      lowerCaseHeaders: true
    }, (err, jsonxls) => {

      async.map(jsonxls, upload_mysql, function (err, results) {
          var number_successful = 0;
          var unsuccessful_entrys = []
          results.forEach(result => {
            if (result[0]) number_successful = number_successful+1;
            else unsuccessful_entrys.push(result[1]);
          });
          callback(number_successful, unsuccessful_entrys);
        });
      });
    };


// Schema that new excel Demo uploads need to follow:
// Image name for corresponding image -- [Probe/id/date]
// File name --> Becomes the Date
// Rest of columns remain.


var upload_mysql = function (row, callback) {
  var values = extract_data_for_probe_sequences_array(row);
  conn.query("INSERT INTO Demo(file_name, link_to_file, slide_name, date," +
    "user, probe, probe_concentration, genotype_a, genotype_b, objective," +
    "optivar, cmount, stages_shown_in_picture, description_of_staining_pattern," +
    "comments, xcoordinate, ycoordinate) VALUES (" + values + ")", (err, res) => {

       if(!err) callback(null, [true, row]);
       else if(err.code == "ER_DUP_ENTRY"){
         winston.error("Duplicate Demo data entered.");
         callback(null, [false, row]);
       }
       else callback(err, [false, row]);
  });
}

function extract_data_for_probe_sequences_array(row) {

// Dynacmically create the image link using other columns or make permanent requirment:
//var link = row['Date'] + row['probe'] + row['name of image'];

var row_to_upload = [];
row_to_upload.push("'" + row['date']+ "'", "'" +row['file name']+"'","'"+row['slide name'] + "'", "'" + row['date']+ "'","'"+ row['user']+"'",
                   "'"+row['probe']+"'", "'" + row['probe concentration'] + "'", "'"+row['genotype a']+"'", "'"+row['genotype b']+"'",
                   "'"+row['objective']+"'","'"+ row['optivar']+"'", "'"+row['c-mount']+"'", "'"+row['stages shown in picture']+"'",
                   "'"+row['description of staining pattern']+"'", "'"+row['comments']+"'", "'"+row['x-coordinate']+"'",
                   "'"+row['y-coordinate']+"'");
                   return row_to_upload;
                 }
