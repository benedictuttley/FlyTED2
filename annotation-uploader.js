//Module to read in annottaions associated with new image uploads.
const xlstojson = require('xls-to-json');
const mysql = require('mysql');
const conn = require('./mysql_setup');
const winston = require('./config/winston'); // Log various actions and store logs in /logs/app.log
const async = require('async');
module.exports.test = function(conn, filename, spreadsheetContent, callback){
console.log("how?");
xlstojson({
      input: filename, // input xls
      output: "output.json", // output json
      lowerCaseHeaders: true
    }, (err, jsonxls) => {
      // Dependent on the type of stylesheet, use the correct sql query.
      console.log("---", spreadsheetContent);
      switch (spreadsheetContent)
      {
        case "image_annotations":
          console.log("33333333333");
          async.mapSeries(jsonxls, upload_mysql, function (err, results) {
            console.log("?????????????");
            var number_successful = 0;
            var unsuccessful_entrys = []
            results.forEach(result => {
              if (result[0]) number_successful = number_successful+1;
              else unsuccessful_entrys.push(result[1]);
            });
            callback(number_successful, unsuccessful_entrys);
          });
          break;

        case "description":
          async.map(jsonxls, upload_mysql_desc, function (err, results) {
            var number_successful = 0;
            var unsuccessful_entrys = []

            results.forEach(result => {
              if (result[0]) number_successful = number_successful+1;
              else unsuccessful_entrys.push(result[1]);
            });

            callback(number_successful, unsuccessful_entrys);

          });
          break;
          default:
        }
      });
    };


// Schema that new excel Demo uploads need to follow:
// Image name for corresponding image -- [Probe/id/date]
// File name --> Becomes the Date
// Rest of columns remain.


var upload_mysql = function (row, cb) {
  console.log("firrred!");
  var values = extract_data_for_probe_sequences_array(row);
  var query = "INSERT INTO Demo(file_name, link_to_file, slide_name, date," +
    "user, probe, probe_concentration, genotype_a, genotype_b, objective," +
    "optivar, cmount, stages_shown_in_picture, description_of_staining_pattern," +
    "comments, xcoordinate, ycoordinate) VALUES (" + values + ")";
  conn.query(query, (err, res) => {
    console.log("pppp");
    if(!err) cb(null, [true, row]);
    else if(err.code == "ER_DUP_ENTRY"){
      winston.error("Duplicate Demo data entered.");
      cb(null, [false, row]);
    }
   else cb(err, [false, row]);
 });
}
//
function extract_data_for_probe_sequences_array(row) {

// Dynacmically create the image link using other columns or make permanent requirment:
//var link = row['Date'] + row['probe'] + row['name of image'];

let row_to_upload = [];
row_to_upload.push("'" + row['date']+ "'", "'" +row['file name']+"'","'"+row['slide name'] + "'", "'" + row['date']+ "'","'"+ row['user']+"'",
                   "'"+row['probe']+"'", "'" + row['probe concentration'] + "'", "'"+row['genotype a']+"'", "'"+row['genotype b']+"'",
                   "'"+row['objective']+"'","'"+ row['optivar']+"'", "'"+row['c-mount']+"'", "'"+row['stages shown in picture']+"'",
                   "'"+row['description of staining pattern']+"'", "'"+row['comments']+"'", "'"+row['x-coordinate']+"'",
                   "'"+row['y-coordinate']+"'");
                   return row_to_upload;
                 }


var upload_mysql_desc = function (row, callback) {
  var values = extract_data_for_description(row);
  var query = "INSERT INTO Probe_Annotations(Gene,Probe_Set,Transcript_ID," +
    "Target_Description,`wt(tin)`,`wt(mip40)(excised)`,`wt(white)`,aly5,comr," +
    "tomb,nht,nxt1,`mip40(sh2)(bsc)`,`mip40(ey)`,`mip40(ey)(bg4)`,wucRNAi,wucRNAi_aly,`mip40(sh2)(comr)`,`mip40(ey)(aly)`,`aly1(27)`,`aly1(18)`,`aly1(18)(btr)`)" +
    "VALUES (" + values + ")";

    conn.query(query, (err, res) => {
      if(!err) callback(null, [true, row]);
      else if(err.code == "ER_DUP_ENTRY"){
        winston.error("Duplicate Demo data entered.");
        callback(null, [false, row]);
      }
     else callback(err, [false, row]);
    });
}


function extract_data_for_description(row){

  let row_to_upload = [];
  row_to_upload.push("'" + row['Gene Symbol']+ "'", "'" +row['Probe Set ID']+"'","'"+row['Transcript ID(Array Design)'] + "'",
                     "'"+row['Target Description'] + "'","'"+row['wt-tin'] + "'",
                     "'" + row['wt-mip40-excised']+ "'","'"+ row['wt-white']+"'","'"+row['aly5']+"'",
                     "'" + row['comr'] + "'", "'"+row['tomb']+"'", "'"+row['nht']+"'", "'"+row['Nxt1']+"'",
                     "'"+ row['mip40(sh2/bsc)']+"'", "'"+row['mip40(ey)']+"'", "'"+row['mip40(ey);bg4']+"'",
                     "'"+row['wucRNAi']+"'", "'"+row['wucRNAi;aly']+"'", "'"+row['mip40(sh2), comr']+"'",
                     "'"+row['mip40(ey);aly']+"'", "'"+row['aly1-27']+"'", "'"+row['aly1-18']+"'", "'"+row['aly1-18-btr']+"'");
                     return row_to_upload;
                   }
