function create_affy_container(gene) {

  // Fetch previous table element:
  var items = document.querySelectorAll("table");
  var lastchild = items[items.length - 1];

  // Create button:
  var chart_affy_button = document.createElement('div');
  var id = "button_affy_" + gene;
  chart_affy_button.setAttribute("id", id);
  chart_affy_button.classList.add('tab2');
  chart_affy_button.innerHTML = ">> VIEW MICROARRAY DATA FOR GENE: " + gene;
  lastchild.parentNode.insertBefore(chart_affy_button, lastchild.nextSibling);
  document.getElementById(id).onclick = function() {
    toggle_microarray_graph(this.id);
  };

  // Create chart container:
  var chart_affy_element = document.createElement('div');
  var id = "chart_affy_" + gene;
  chart_affy_element.setAttribute("id", id);
  chart_affy_button.appendChild(chart_affy_element);
}


function create_tissue_container(gene) {
  var items = document.querySelectorAll("div"); // Fetch previous table element.
  var lastchild = items[items.length - 1];

  var chart_button = document.createElement('div'); // Create button.
  var id = gene;
  chart_button.setAttribute("id", id);
  chart_button.classList.add('tab');
  chart_button.innerHTML = ">> EXPRESSION BY TISSUE GRAPH FOR GENE: " + gene;
  lastchild.parentNode.insertBefore(chart_button, lastchild.nextSibling);
  document.getElementById(id).onclick = function() {
    toggle_tissue_graph(this.id);
  };

  var chart_element = document.createElement('div'); // Create chart container.
  var id = "chart_" + gene;
  chart_element.setAttribute("id", id);
  chart_button.appendChild(chart_element);
}
