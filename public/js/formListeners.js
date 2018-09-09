/*When user supplied gene name is submitted for quick change,
 fetch each dataset and display results in the editor frame:
 */
var current_gene;
$('#gene_for_change').on('submit', function(e) {
  console.log($('#quick_edit_gene'));
  current_gene = $('#quick_edit_gene').val();
  e.preventDefault();
  document.getElementById('pls').innerHTML = "";
  document.getElementById('myDivGeneInput').innerHTML = "";

  // Check form input is not empty, display error and prevent submission:
  if ($('#quick_edit_gene').val() == "") {
    $('#myDivGeneInput').append("<ul id='errorGeneInput' class='errors'></ul>");
    $("#errorGeneInput").append("<li>Gene field is empty</li>");

    // Clear all input fields on bad gene input:
    var input_fileds = document.getElementsByTagName("input");
    for (var i = 0; i < input_fileds.length; i++) {
      if (input_fileds[i].type == "text") {
        input_fileds[i].value = "";
      }
    }
  }

  // Fetch the microarray data and display in the microarray editor frame:
  else {
    $.ajax({
      type: 'post',
      url: '/quick_change_microarray_fetch',
      data: {
        gene: $('#quick_edit_gene').val()
      },
      success: function(data) {
        console.log("***");
        console.log(data[0]);
        console.log("***");
        for (var col in data[0]) {
          if (data[0].hasOwnProperty(col)) {
            // Add data value to corresponding column input field,
            // if column value is empty leave input field empty: */
            if (document.getElementById(col) != null) {
              document.getElementById(col).value = data[0][col];
            }
          }
        }
      }
    });

    // Fetch the probe sequence data and display in the editor frame:
    console.log($('#quick_edit_gene'));
    $.ajax({
      type: 'post',
      url: '/quick_change_probe_fetch',
      data: {
        gene: $('#quick_edit_gene').val()
      },
      success: function(data) {
        for (var col in data) {
          if (data.hasOwnProperty(col)) {
            // Add data value to corresponding column input field,
            // if column value is empty leave input field empty: */
            if (document.getElementById(col) != null) {
              document.getElementById(col).value = data[col];
            }
          }
        }
      }

    });


    // Fetch the image metadata and display in the editor frame:
    $.ajax({
      type: 'post',
      url: '/quick_change_image_fetch',
      data: {
        gene: $('#quick_edit_gene').val()
      },
      success: function(data) {
        // Reset image metadata view:
        $('.your-class').replaceWith("<div class='your-class' id='pls'></div>");

        var counter = 0;
        document.getElementById('pls').innerHTML = "";
        document.getElementById('myDivseven').innerHTML = ""; // clear quick image notifier.
        document.getElementById('myDivfive').innerHTML = ""; // Clear microarray notifier.
        document.getElementById('myDivsix').innerHTML = ""; // Clear probe notifier.

        displayAdditionContainer("image_container_tab", true);
        displayAdditionContainer("microarray_container_tab", true);
        displayAdditionContainer("probe_container_tab", true);

        data.forEach(set => {

          // Create form
          var idforthis = "test" + counter;
          var f = document.createElement("form");
          f.setAttribute('id', idforthis);
          f.setAttribute('class', "quick_change_form")
          f.enctype = 'multipart/form-data';
          f.method = 'post';
          var il = document.createElement('div');
          il.setAttribute('class', 'form-group row');

          // Create delete button:
          let del_button = document.createElement("button");
          del_button.innerHTML = "Delete Dataset";
          del_button.setAttribute('type', 'button');
          del_button.setAttribute('class', 'btn');
          del_button.id = set.link_to_file;
          del_button.setAttribute('onclick',
            'delImgAndData(this.id)');
          f.appendChild(del_button);

          for (var col in set) {
            if (set.hasOwnProperty(col)) {
              console.log(col + " -> " + set[col]);
              var p = document.createElement("div");
              p.setAttribute('class', "col-xs-2");
              var i = document.createElement("input"); //input element, text
              i.setAttribute('class', "form-control");
              i.setAttribute('type', "text");
              i.setAttribute('name', col);
              i.value = set[col];
              var l = document.createElement("label"); //input element, Submit button
              l.innerHTML = col;
              p.appendChild(l);
              p.appendChild(i);
              il.appendChild(p);
            }
          }

          var t = document.createElement("input"); //Hidden for orignial [1]
          t.setAttribute('type', "text");
          t.setAttribute('name', 'hidden_original_link');
          t.value = set.link_to_file;
          t.style.display = "none";
          il.appendChild(t);

          var n = document.createElement("input"); //Hidden for orignial [2]
          n.setAttribute('type', "text");
          n.setAttribute('name', 'hidden_original_probe');
          n.value = set.probe;
          n.style.display = "none";
          il.appendChild(n);
          f.append(il);

          // Create div to contain a image/metadata unit:
          let test_container = document.createElement('div');

          // Display image:
          let img = document.createElement("img");
          img.src = set.link_to_file;
          img.style.position = "relative";
          img.style.width = "50%";
          img.style.height = "50%";
          img.setAttribute('id', 'ok_nice');
          test_container.append(img)

          f.appendChild(test_container);
          document.getElementById('pls').appendChild(f);
          f.appendChild(document.createElement('br'));

          // Create form submission button:
          let form_sumbit_button = document.createElement("button");
          form_sumbit_button.setAttribute("type", "submit");
          form_sumbit_button.setAttribute('class', 'btn');
          form_sumbit_button.innerHTML = "Update Dataset";
          f.append(form_sumbit_button);

          // Create and display custom file submit input:
          var b = set.link_to_file;
          var input = (`<input type="file" name="yo_${b}" id="yo_${b}" class="inputfile">`);
          var label = (
            `<label for="yo_${b}"><svg xmlns="http://www.w3.org/` +
            `2000/svg" width="20" height="17" viewBox="0 0 20 17"><path d="M10 0l-5` +
            `.2 4.9h3.3v5.1h3.8v-5.1h3.3l-5.2-4.9zm9.3 11.5l-3.2-2.1h-2l3.4 2.6h-3.` +
            `5c-.1 0-.2.1-.2.1l-.8 2.3h-6l-.8-2.2c-.1-.1-.1-.2-.2-.2h-3.6l3.4-2.6h-` +
            `2l-3.2 2.1c-.4.3-.7 1-.6 1.5l.6 3.1c.1.5.7.9 1.2.9h16.3c.6 0 1.1-.4 1.` +
            `3-.9l.6-3.1c.1-.5-.2-1.2-.7-1.5z"/></svg> <span>Select Image Replaceme` +
            `nt&hellip;</span></label>`);

          var t = document.getElementById('content');

          var temp = document.createElement('div');
          temp.innerHTML = input;

          var temptwo = document.createElement('div');
          temptwo.innerHTML = label;

          var p = document.createElement('div');
          p.id = "p";
          p.appendChild(document.createElement('br'));

          test_container.appendChild(p);
          var t = document.createElement('div');
          var e = document.createElement('div');
          t.appendChild(temp.firstChild);
          t.appendChild(temptwo.firstChild);
          p.appendChild(t);
          p.appendChild(e);
          customFileInput()

          // Fetch the form that is being submitted out of all forms in the slder:
          $('#' + idforthis).on('submit', function(e) {
            document.getElementById('myDivseven').innerHTML = "";
            e.preventDefault();

            // Get form and create an FormData object:
            let form = $('#' + idforthis)[0];
            var data = new FormData(form);
            document.getElementById('ok_nice').src = "";

            // Submit changes to image metadata or the image file itself:
            $.ajax({
              type: "POST",
              enctype: 'multipart/form-data',
              url: "/quick_change_img_edit",
              data: data, // pass the form data object as data to be uploaded
              processData: false,
              contentType: false,
              cache: false,
              timeout: 600000,
              success: function(response) {

                // Used to generate unique image file path that will ensure image is reloaded
                let d = new Date();
                document.getElementById('ok_nice').src = form.link_to_file.value + "?ver="; + d.getTime();
                // On upload completion server side, display confirmation message:
                if (response.success != undefined && response.success.length > 0) {
                  $('#myDivseven').append(
                    "<ul id='successListSeven' class='success'></ul>"
                  );
                  for (let i = 0; i < response.success.length; i++) {
                    $("#successListSeven").append("<li>" +
                      response.success[i] + "</li>");
                  }
                }
              }
            });
          });
          counter++;
        });

        // Generate slider with each slide being an editor for one image metadata unit:
        $('.your-class').slick({
          infinite: false,
          swipe: false,
          prevArrow: $('#prev_slider'),
          nextArrow: $('#next_slider')
        });
      }
    });
  }

});

// Submit changes to microarray data:
$('#microarray_form').on('submit', function(e) {
  document.getElementById('myDivfive').innerHTML = "";
  e.preventDefault();
  document.getElementById('test_Gene').value = current_gene;

  $.ajax({
    type: 'post',
    url: '/quick_change_microarray_edit',
    data: $(this).serialize(),
    success: function(response) {
      // On upload completion server side, display confirmation message:
      if (response.success != undefined && response.success.length > 0) {

        $('#myDivfive').append(
          "<ul id='successListFive' class='success'></ul>");

        for (let i = 0; i < response.success.length; i++) {
          $("#successListFive").append("<li>" + response.success[i] +
            "</li>");
        }
      }
    }
  });
});


// Submit new image and metadata unit:
$('#image_addition_form').on('submit', function(e) {
  e.preventDefault();
  var form = $('#image_addition_form')[0];
  var data = new FormData(form);
  document.getElementById('image_addition_notifier').innerHTML = "";
  $.ajax({
    type: 'post',
    enctype: 'multipart/form-data',
    url: '/QuickImageAddition',
    data: data,
    processData: false,
    contentType: false,
    cache: false,
    timeout: 600000,
    success: function(response) {

      // On upload completion server side, display error message, if error recieved:
      if (response.error != undefined && response.error.length > 0) {
        $('#image_addition_notifier').append(
          "<ul id='errorListAddition' class='errors'></ul>");
        for (let i = 0; i < response.error.length; i++) {
          $("#errorListAddition").append("<li>" + response.error[i] + "</li>");
        }
      }
      // On upload completion server side, display confirmation message, if success recieved:
      if (response.success != undefined && response.success.length > 0) {
        $('#image_addition_notifier').append(
          "<ul id='successListAddition' class='success'></ul>");
        for (let i = 0; i < response.success.length; i++) {
          $("#successListAddition").append("<li>" + response.success[i] +
            "</li>");
        }
      }
    }
  });
});

// Submit changes to probe data:
$('#probe_form').on('submit', function(e) {
  e.preventDefault();
  document.getElementById('myDivsix').innerHTML = "";

  // Show loaing div when performing probe Blast as it takes some time:
  $('#myDivsix').append(
    "<ul id='loadingListsix' class='loading'></ul>"
  );
  $("#loadingListsix").append("<li>Please Wait...</li>");

  $.ajax({
    type: 'post',
    url: '/quick_change_probe_edit',
    data: $(this).serialize(),
    success: function(response) {
      document.getElementById('myDivsix').innerHTML = "";


      // On upload completion server side, display error message, if error recieved:
      if (response.error != undefined && response.error.length > 0) {
        $('#myDivsix').append(
          "<ul id='errorListSix' class='errors'></ul>");
        for (let i = 0; i < response.error.length; i++) {
          $("#errorListSix").append("<li>" + response.error[i] + "</li>");
        }
      }
      // On upload completion server side, display confirmation message, if success recieved:
      if (response.success != undefined && response.success.length > 0) {
        $('#myDivsix').append(
          "<ul id='successListSix' class='success'></ul>");
        for (let i = 0; i < response.success.length; i++) {
          $("#successListSix").append("<li>" + response.success[i] +
            "</li>");
        }
      }
    }
  });
});


function customFileInput() {
  ///////////////////////////////////////////////////////////////////////
  /*
    By Osvaldas Valutis, www.osvaldas.info
    Available for use under the MIT License
  */

  'use strict';

  ;
  (function(document, window, index) {
    var inputs = document.querySelectorAll('.inputfile');
    Array.prototype.forEach.call(inputs, function(input) {
      var label = input.nextElementSibling,
        labelVal = label.innerHTML;

      input.addEventListener('change', function(e) {
        var fileName = '';
        if (this.files && this.files.length > 1)
          fileName = (this.getAttribute(
            'data-multiple-caption') || '').replace(
            '{count}', this.files.length);
        else
          fileName = e.target.value.split('\\').pop();

        if (fileName)
          label.querySelector('span').innerHTML =
          fileName;
        else
          label.innerHTML = labelVal;
      });

      // Firefox bug fix
      input.addEventListener('focus', function() {
        input.classList.add('has-focus');
      });
      input.addEventListener('blur', function() {
        input.classList.remove('has-focus');
      });
    });
  }(document, window, 0));
}
