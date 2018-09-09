// Module to read in annotations associated with new image uploads.
const xlstojson = require('xls-to-json');
const mysql = require('mysql');
const conn = require('./mysql_setup');
const winston = require('./config/winston'); // Log various actions and store logs in /logs/app.log
const async = require('async');

module.exports.test = function(conn, filename, spreadsheetContent, callback) {
  xlstojson({
    input: filename, // input xls
    output: "output.json", // output json
    lowerCaseHeaders: true
  }, (err, jsonxls) => {

    // Dependent on the type of stylesheet, use the correct sql query.
    switch (spreadsheetContent) {
      case "image_annotations":

        async.mapSeries(jsonxls, upload_mysql, function(err, results) {

          var number_successful = 0;
          var unsuccessful_entrys = []
          results.forEach(result => {
            if (result[0]) number_successful = number_successful + 1;
            else unsuccessful_entrys.push(result[1]);
          });
          callback(number_successful, unsuccessful_entrys);
        });
        break;

      case "description":
        async.map(jsonxls, upload_mysql_desc, function(err, results) {
          var number_successful = 0;
          var unsuccessful_entrys = []

          results.forEach(result => {
            if (result[0]) number_successful = number_successful + 1;
            else unsuccessful_entrys.push(result[1]);
          });

          callback(number_successful, unsuccessful_entrys);

        });
        break;
      default:
    }
  });
};


var upload_mysql = function(row, cb) {
  var values = extract_data_for_probe_sequences_array(row);
  let query = "INSERT INTO Demo(file_name, link_to_file, slide_name, date," +
    "user, probe, probe_concentration, genotype_a, genotype_b, objective," +
    "optivar, cmount, stages_shown_in_picture, description_of_staining_pattern," +
    "comments, xcoordinate, ycoordinate) VALUES (" + values + ")";
  conn.query(query, (err, res) => {
    if (!err) cb(null, [true, row]);
    else if (err.code == "ER_DUP_ENTRY") {
      winston.error("Duplicate Demo data entered.");
      cb(null, [false, row]);
    } else cb(err, [false, row]);
  });
}

var upload_mysql_desc = function(row, callback) {
  var values = extract_data_for_description(row);
  var query = "INSERT INTO Probe_Annotations(Gene,Probe_Set,Transcript_ID," +
  "Target_Description,`wt(tin)`,`wt(mip40)(excised)`,`wt(white)`,aly5,comr," +
  "tomb,nht,nxt1,`mip40(sh2)(bsc)`,`mip40(ey)`,`mip40(ey)(bg4)`,wucRNAi," +
  "wucRNAi_aly,`mip40(sh2)(comr)`,`mip40(ey)(aly)`,`aly1(27)`,`aly1(18)`," +
  "`aly1(18)(btr)`)" + "VALUES (" + values + ")";

  conn.query(query, (err, res) => {
    if (!err) callback(null, [true, row]);
    else if (err.code == "ER_DUP_ENTRY") {
      winston.error("Duplicate Demo data entered.");
      callback(null, [false, row]);
    } else callback(err, [false, row]);
  });
}

function extract_data_for_probe_sequences_array(row) {

  let row_to_upload = [];
  row_to_upload.push(conn.escape(row['date']), conn.escape(row['file name']), conn.escape(row['slide name']),
  conn.escape(row['date']), conn.escape(row['user']), conn.escape(row['probe']), conn.escape(row['probe concentration']),
  conn.escape(row['genotype a']), conn.escape(row['genotype b']),
  conn.escape(row['objective']), row['optivar'], row['c-mount'], row['stages shown in picture'],
  conn.escape(row['description of staining pattern']), conn.escape(row['comments']), conn.escape(row['x-coordinate']),
  conn.escape(row['y-coordinate']));

  return row_to_upload;
}

function extract_data_for_description(row) {

  let row_to_upload = [];
  row_to_upload.push(conn.escape(row['Gene Symbol']),conn.escape(row['Probe Set ID']),conn.escape(row['Transcript ID(Array Design)']),
  conn.escape(row['Target Description']),conn.escape(row['wt-tin']),conn.escape(row['wt-mip40-excised']),conn.escape(row['wt-white']),
  conn.escape(row['aly5']),conn.escape(row['comr']),conn.escape(row['tomb']),conn.escape(row['nht']),conn.escape(row['Nxt1']),
  conn.escape(row['mip40(sh2/bsc)']),conn.escape(row['mip40(ey)']),conn.escape(row['mip40(ey);bg4']),conn.escape(row['wucRNAi']),
  conn.escape(row['wucRNAi;aly']),conn.escape(row['mip40(sh2), comr']),conn.escape(row['mip40(ey);aly']),conn.escape(row['aly1-27']),
  conn.escape(row['aly1-18']),conn.escape(row['aly1-18-btr']));

  return row_to_upload;
}


module.exports.quickAdditionImageUpload = function(values, cb) {
  let query = "INSERT INTO Demo(file_name, link_to_file, slide_name, date," +
    "user, probe, probe_concentration, genotype_a, genotype_b, objective," +
    "optivar, cmount, stages_shown_in_picture, description_of_staining_pattern," +
    "comments, xcoordinate, ycoordinate) VALUES (" + values + ")";
  conn.query(query, (err, res) => {
    if (!err) {
      cb(null);
    } else if (err.code == "ER_DUP_ENTRY") {
      winston.error("Duplicate Demo data entered.");
      cb(err);
    } else {
      winston.error(`${err.status || 500} - ${err.message}`);
      cb(err);
    }
  });

}
