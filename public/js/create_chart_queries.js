// --- SET OPTIONS AND DRAW CHART FOR AFFY DATA START ---

function query_for_affy_chart(affy_data, gene) {
  var data = google.visualization.arrayToDataTable(affy_data);
  var options = { // Set the options for the affy chart.
    legend: 'bottom',
    width: '100%',
    height: 800,
    title: 'Microarray Results For Whole Testes',
    colors: ['#4286f4'],
    hAxis: {
      title: 'Mutant',
      slantedText: true,
      slantedTextAngle: 90,
    },
    chartArea: {
      width: "82%",
      top: 28,
      height: '60%'
    },
    vAxis: {
      title: 'Expression Count',
      format: 'decimal',
      slantedText: true,
      slantedTextAngle: 90,
      sliceVisibilityThreshold: 0
    },
  };

  let id = "chart_affy_" + gene;
  var container = document.getElementById(id)
  var chart_affy = new google.visualization.ColumnChart(container);
  container.style.display = 'block';
  google.visualization.events.addListener(chart_affy, 'ready', () => {
    container.style.display = 'none'
  });

  affy_mine.push(gene, data, options);
  chart_affy.draw(data, options);
}
// --- SET OPTIONS AND DRAW CHART FOR AFFY DATA END ---



// --- SET OPTIONS AND DRAW CHART FOR EXPRESSION DATA START ---

function testing(data_array, gene) {
  var data = google.visualization.arrayToDataTable(data_array);
  var options = { // Set the options for the tissue chart.
    legend: 'bottom',
    width: '100%',
    height: 800,
    title: 'Gene Expression in Different Drosopilla Tissues - Data From FlyAtlas',
    colors: ['#d4374a'],
    hAxis: {
      title: 'Tissue Type',
      slantedText: true,
      slantedTextAngle: 90,
    },
    chartArea: {
      width: "82%",
      top: 28,
      height: '60%'
    },
    vAxis: {
      title: 'Expression (mRNA Signal SEM)',
      format: 'short',
      slantedText: true,
      slantedTextAngle: 90,
      sliceVisibilityThreshold: 0
    },
  };
  let id = "chart_" + gene;
  var container = document.getElementById(id)
  var chart = new google.visualization.ColumnChart(container);
  container.style.display = 'block';
  google.visualization.events.addListener(chart, 'ready', () => {
    container.style.display = 'none'
  });

  mine.push(gene, data, options);
  chart.draw(data, options);
}
// --- SET OPTIONS AND DRAW CHART FOR EXPRESSION DATA END ---
