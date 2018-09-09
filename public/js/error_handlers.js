function Handle_No_Tissue_Data(probe) {
  let tissue_title = document.getElementById("expression_title_" + probe);
  tissue_title.innerHTML = ">>> NO EXPRESSION DATA AVAILABLE!";
  document.getElementById(probe + "_table").remove(); // Remove tissue expression table.
  document.getElementById("expression_form_" + probe).remove(); // Remove tissue expression form
}

function Handle_No_Microarray_Data(probe) {
  let affy_title = document.getElementById("microarray_" + probe);
  affy_title.innerHTML = ">>> NO MICROARRAY DATA AVAILABLE!";
  document.getElementById("annotation_for_" + probe).remove();
}

function Handle_No_Image_Data(probe) {
  let image_title = document.getElementById(probe + "_image_title");
  image_title.innerHTML = ">>> NO IMAGE DATA AVAILABLE!";
}

function Check_For_Data_Components(microarray_data, tissue_data, transcript_data, flyted_data) {
  console.log(flyted_data);
  if (microarray_data === undefined || microarray_data.length == 0) microarray_data_present = false; //CHECK IF MICROARRAY DATA EXIST FOR THIS GENE.
  if (tissue_data === undefined || tissue_data.length == 0) tissue_expression_present = false; //CHECK IF TISSUE EXPRESSION DATA EXIST FOR THIS GENE.
  if (transcript_data === undefined || transcript_data.length == 0) transcipt_expression_present = false; //CHECK IF TRANSCRIPT DATA EXIST FOR THIS GENE.
  if (flyted_data === undefined || flyted_data.length == 0) image_data_present = false; //CHECK IF IMAGES AND ASSOCIATED METADATA EXIST FOR THIS GENE.
}
