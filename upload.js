// Imports used in module
const app = require('./start');
const fs = require('fs');
const mkdirp = require('mkdirp');
const multer = require('multer');
const winston = require('./config/winston'); // Log various actions and store logs in /logs/app.log
const storage = multer.diskStorage({
  destination: function (req, file, cb) {
    cb(null, 'multerTemp')
  },
  filename: function (req, file, cb) {
    cb(null, file.originalname)
  }
});
const testmulter = multer({ storage: storage })

module.exports.read = function (req, res, callback) {
  var counter = 0;
  var files = req.files;
  var num_uploaded = 0;

  files.forEach((file) => {
    // Assign the image descriptors
    var descriptors = file.originalname.split("_");
   console.log(descriptors);
    var image_gene = descriptors[0];
    var image_date = descriptors[2].split('.').join('-').replace("-bmp", "");
    var parent_directory = ("/var/www/html/FlyTED2/public/img/" + image_date + "/" + image_gene);  // Root directory for all uploaded images.
    var complete_file_location = parent_directory + "/" + file.originalname;  // Target directory for file.
    console.log(complete_file_location);
    mkdirp(parent_directory, err => { // Make correct directory based on file name

    if (err) winston.error(`${err.status || 500} - ${err.message} - ${req.originalUrl} - ${req.method} - ${req.ip}`);

    fs.access(complete_file_location, fs.constants.R_OK, (err) => {
      if(err){
        copyImage(file, complete_file_location, err => {
          counter = counter + 1;
          if (err) winston.error(`${err.status || 500} - ${err.message} - ${req.originalUrl} - ${req.method} - ${req.ip}`);
          else {
            req.flash('success', file.originalname + " uploaded successfully");
            num_uploaded = num_uploaded + 1;
            let completed = checkForUploadCompletion(counter, res, files.length);
            if (completed) callback(num_uploaded);
          }  // If file did not previously exist, add success notification.

        // removeOld(file.path, err => {
        //   if (err) winston.error(`${err.status || 500} - ${err.message} - ${req.originalUrl} - ${req.method} - ${req.ip}`);
        //   let completed = checkForUploadCompletion(counter, res, files.length);
        //   if (completed) callback(num_uploaded);
        // });
      });
    }

    else {
        counter = counter + 1;
        winston.warn(`File ${file.originalname} already exists, preventing overwrite`);
        req.flash('error',file.originalname + " upload failed, file already exists"); // If file already exists add error notification.
        let completed = checkForUploadCompletion(counter, res, files.length);
        if (completed) callback(num_uploaded);
      }
    });
  });
  });
}


// Cleaner Method:
function removeOld(file_path, callback){
  fs.unlink(file_path, (err) => { // When files uploaded, remove them from temp storage.
  if (err) callback(err);
  else callback(null);
  });
}

// Return notifications on completion
function checkForUploadCompletion(counter, res, final_length){
  if(counter >= final_length){  // When all async calls have completed, return notification array containing success and error objects for view display.
    return true;
  } else {
    return false;
  }
}

// Method to permorm file copy to correct directory
function copyImage(file, complete_file_location, callback) {
  fs.rename(file.path, complete_file_location, err => {
    console.log("image should be copied");
  if(err) callback(err);
  else callback(null);
  });
}


module.exports.uploadexcel = function (req, res, callback) {
  callback(req.files[0].path);
}
