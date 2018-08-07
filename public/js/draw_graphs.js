// --- FUNCTION TO DRAW CHARTS USING EXPRESSION AND MICROARRAY DATA START ---


function drawChart() { // For each gene in the original query, draw the respective tissue expression and affy value graph.
  if (tissue_expression_present) {
    for (let i = 1; i < data_array.length; i = i + 2) {
      testing(data_array[i], data_array[i - 1]);
    }
  }
  if (microarray_data_present) {
    for (let i = 1; i < affy_array.length; i = i + 2) {
      query_for_affy_chart(affy_array[i], affy_array[i - 1]);
    }
  }
}
// --- FUNCTION TO DRAW CHARTS USING EXPRESSION AND MICROARRAY DATA END--


// --- FUNCTIONS TO TOGGLE DISPLAY OF CHARTS START ---

function toggle_tissue_graph(element) {

  var y = element.replace("chart_", "");
  var z = "chart_" + y;
  var x = document.getElementById(z);
  if (x.style.display === "none") {
    x.style.display = "block";
  } else {
    x.style.display = "none";
  }
}

function toggle_microarray_graph(element) { // Toogle the display of the affy chart when the tab above the chart is clicked.

  var y = element.replace("button_affy_", "chart_affy_");
  var x = document.getElementById(y);
  if (x.style.display === "none") {
    x.style.display = "block";
  } else {
    x.style.display = "none";
  }
}

// --- FUNCTIONS TO TOGGLE DISPLAY OF CHARTS END ---

// --- FUNCTIONS TO REDRAW CHART ON WINDOW DIMENSION CHANGE START ---

function resize_graphs_on_screen_change() { // When window is resized, redraw every expression chart to fit the new window dimensions.
  for (let i = 0; i < mine.length; i += 3) {
    let hidden_before = false;
    let thisChart = document.getElementById("chart_" + mine[i]);
    if (thisChart.style.display == "none") hidden_before = true;

    let chart = new google.visualization.ColumnChart(thisChart);
    thisChart.style.display = 'block';
    google.visualization.events.addListener(chart, 'ready', () => {
      if (hidden_before) thisChart.style.display = 'none'
    });
    chart.draw(mine[i + 1], mine[i + 2]);
  }
  // When window is resized, redraw every affy chart to fit the new window dimensions.
  for (let i = 0; i < affy_mine.length; i += 3) {
    let hidden_before = false;
    let thisChart = document.getElementById("chart_affy_" + mine[i]);
    if (thisChart.style.display == "none") hidden_before = true;

    let chart = new google.visualization.ColumnChart(thisChart);
    thisChart.style.display = 'block';
    google.visualization.events.addListener(chart, 'ready', () => {
      if (hidden_before) thisChart.style.display = 'none'
    });
    chart.draw(affy_mine[i + 1], affy_mine[i + 2]);
  }
}

// --- FUNCTIONS TO REDRAW CHART ON WINDOW DIMENSION CHANGE END ---
