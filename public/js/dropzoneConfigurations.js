// Configure the Dropzone.
Dropzone.options.myAwesomeDropzone = {

  //autoProcessQueue: false,
  uploadMultiple: true,
  parallelUploads: 100,

  init: function() {
    var myDropzone = this;

    this.on("drop", function(e) {
      document.getElementById('myDiv').innerHTML = "";
      e.preventDefault();
      e.stopPropagation();
      myDropzone.processQueue();
    });


    // When server has finished proce    <button type="submit" class="image_submit">Submit</button>ssing, a success message is recieved and notifications are displayed
    this.on("successmultiple", (files, response) => {
      document.getElementById('my-awesome-dropzone').innerHTML = "";
      document.getElementById('my-awesome-dropzone').innerHTML =
        '<div class="dz-message" data-dz-message="">' +
        '<span style="font-weight: bold">Submit Expression Images Here: Read the Descripton</span></div>';
      document.getElementById('myDiv').innerHTML = "";

      if (response.error != undefined && response.error.length > 0) {
        $('#myDiv').append("<ul id='errorList' class='errors'></ul>");
        for (let i = 0; i < response.error.length; i++) {
          $("#errorList").append("<li>" + response.error[i] + "</li>");
        }
      }

      if (response.success != undefined && response.success.length > 0) {
        $('#myDiv').append("<ul id='successList' class='success'></ul>");
        for (let i = 0; i < response.success.length; i++) {
          $("#successList").append("<li>" + response.success[i] +
            "</li>");
        }
      }
    });

    this.on("errormultiple", function(files, response) {
      document.getElementById('my-awesome-dropzone').innerHTML = "";
      document.getElementById('my-awesome-dropzone').innerHTML =
        '<span style="font-weight: bold">Submit Expression Images Here: Read the Descripton</span></div>';
      document.getElementById('myDiv').innerHTML = "";
      $('#myDiv').append("<ul id='errorList' class='errors'></ul>");
      $("#errorList").append("<li>" +
        "Error connecting to server, it may be down" + "</li>");


      // Gets triggered when there was an error sending the files.
      // Maybe show form again, and notify user of error
    });
  }
}



Dropzone.options.myAwesomeDropzoneTwo = {

  //autoProcessQueue: false,
  uploadMultiple: true,
  parallelUploads: 100,

  init: function() {
    var myDropzone = this;

    this.on("drop", function(e) {
      document.getElementById('myDivtwo').innerHTML = "";
      e.preventDefault();
      e.stopPropagation();
      myDropzone.processQueue();
    });


    // When server has finished proce    <button type="submit" class="image_submit">Submit</button>ssing, a success message is recieved and notifications are displayed
    this.on("successmultiple", (files, response) => {

      document.getElementById('my-awesome-dropzone-two').innerHTML = "";
      document.getElementById('my-awesome-dropzone-two').innerHTML =
        '<div class="dz-message" data-dz-message="">' +
        '<span style="font-weight: bold">Submit Expression Images Here: Read the Descripton</span></div>';
      document.getElementById('myDivtwo').innerHTML = "";

      if (response.error != undefined && response.error.length > 0) {
        $('#myDivtwo').append("<ul id='two_errorList' class='errors'></ul>");
        for (let i = 0; i < response.error.length; i++) {
          $("#two_errorList").append("<li>" + response.error[i] + "</li>");
        }
      }

      if (response.success != undefined && response.success.length > 0) {
        $('#myDivtwo').append(
          "<ul id='two_successList' class='success'></ul>");
        for (let i = 0; i < response.success.length; i++) {
          $("#two_successList").append("<li>" + response.success[i] +
            "</li>");
        }
      }
    });

    this.on("errormultiple", function(files, response) {
      document.getElementById('my-awesome-dropzone-two').innerHTML = "";
      document.getElementById('my-awesome-dropzone-two').innerHTML =
        '<span style="font-weight: bold">Submit Expression Images Here: Read the Descripton</span></div>';
      document.getElementById('myDivtwo').innerHTML = "";
      $('#myDivtwo').append("<ul id='two_errorList' class='errors'></ul>");
      $("#two_errorList").append("<li>" +
        "Error connecting to server, it may be down" + "</li>");


      // Gets triggered when there was an error sending the files.
      // Maybe show form again, and notify user of error
    });
  }
}


Dropzone.options.myAwesomeDropzoneThree = {

  //autoProcessQueue: false,
  uploadMultiple: true,
  parallelUploads: 100,

  init: function() {
    var myDropzone = this;

    this.on("drop", function(e) {
      document.getElementById('myDivthree').innerHTML = "";
      e.preventDefault();
      e.stopPropagation();
      myDropzone.processQueue();
    });


    // When server has finished proce    <button type="submit" class="image_submit">Submit</button>ssing, a success message is recieved and notifications are displayed
    this.on("successmultiple", (files, response) => {

      document.getElementById('my-awesome-dropzone-three').innerHTML =
        "";
      document.getElementById('my-awesome-dropzone-three').innerHTML =
        '<div class="dz-message" data-dz-message="">' +
        '<span style="font-weight: bold">Submit Expression Images Here: Read the Descripton</span></div>';
      document.getElementById('myDivthree').innerHTML = "";

      if (response.error != undefined && response.error.length > 0) {
        $('#myDivthree').append(
          "<ul id='three_errorList' class='errors'></ul>");
        for (let i = 0; i < response.error.length; i++) {
          $("#three_errorList").append("<li>" + response.error[i] + "</li>");
        }
      }

      if (response.success != undefined && response.success.length > 0) {
        $('#myDivthree').append(
          "<ul id='three_successList' class='success'></ul>");
        for (let i = 0; i < response.success.length; i++) {
          $("#three_successList").append("<li>" + response.success[i] +
            "</li>");
        }
      }
    });

    this.on("errormultiple", function(files, response) {
      document.getElementById('my-awesome-dropzone-three').innerHTML =
        "";
      document.getElementById('my-awesome-dropzone-three').innerHTML =
        '<span style="font-weight: bold">Submit Description spreadsheets here.</span></div>';
      document.getElementById('#myDivthree').innerHTML = "";
      $('#myDivthree').append("<ul id='three_errorList' class='errors'></ul>");
      $("#three_errorList").append("<li>" +
        "Error connecting to server, it may be down" + "</li>");
    });
  }
}


Dropzone.options.myAwesomeDropzoneFour = {

  //autoProcessQueue: false,
  uploadMultiple: true,
  parallelUploads: 100,

  init: function() {
    var myDropzone = this;

    this.on("drop", function(e) {
      document.getElementById('myDivfour').innerHTML = "";
      e.preventDefault();
      e.stopPropagation();
      myDropzone.processQueue();
    });


    // When server has finished proce    <button type="submit" class="image_submit">Submit</button>ssing, a success message is recieved and notifications are displayed
    this.on("successmultiple", (files, response) => {

      document.getElementById('my-awesome-dropzone-four').innerHTML =
        "";
      document.getElementById('my-awesome-dropzone-four').innerHTML =
        '<div class="dz-message" data-dz-message="">' +
        '<span style="font-weight: bold">Submit Expression Images Here: Read the Descripton</span></div>';
      document.getElementById('myDivfour').innerHTML = "";


      if (response.error != undefined && response.error.length > 0) {
        $('#myDivfour').append(
          "<ul id='four_errorList' class='errors'></ul>");
        for (let i = 0; i < response.error.length; i++) {
          $("#four_errorList").append("<li>" + response.error[i] + "</li>");
        }
      }

      if (response.success != undefined && response.success.length > 0) {
        $('#myDivfour').append(
          "<ul id='five_successList' class='success'></ul>");
        for (let i = 0; i < response.success.length; i++) {
          $("#five_successList").append("<li>" + response.success[i] +
            "</li>");
        }
      }
    });

    this.on("errormultiple", function(files, response) {
      document.getElementById('my-awesome-dropzone-four').innerHTML =
        "";
      document.getElementById('my-awesome-dropzone-four').innerHTML =
        '<span style="font-weight: bold">Submit Description spreadsheets here.</span></div>';
      document.getElementById('#myDivfour').innerHTML = "";
      $('#myDivfour').append("<ul id='four_errorList' class='errors'></ul>");
      $("#four_errorList").append("<li>" +
        "Error connecting to server, it may be down" + "</li>");
    });
  }
}
