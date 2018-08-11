// --- ADD ROW FILTERING TO BY TISSUE TO TISSUE EXPRESSION TABLE START ---
function SearchTable(id) {
  var table_id, input, filter, table, tr, td, i;
  input = document.getElementById(id);
  filter = input.value.toUpperCase();

  // Construct the table ID:
  table_id = id.replace('input', 'table');
  table = document.getElementById(table_id);
  tr = table.getElementsByTagName("tr");
  for (i = 0; i < tr.length; i++) {
    td = tr[i].getElementsByTagName("td")[8];
    if (td) {
      if (td.innerHTML.toUpperCase().indexOf(filter) > -1) {
        tr[i].style.display = "";
      } else {
        tr[i].style.display = "none";
      }
    }
  }
}
// --- ADD ROW FILTERING TO BY TISSUE TO TISSUE EXPRESSION TABLE END ---

// --- SORT THE TABLE IN ASCENDING/DESCENDING ORDER, CREDIT: W3C. START ---
function sortTable(gene, order) {
  var table, rows, switching, i, x, y, shouldSwitch;
  table = document.getElementById((gene + "_table"));
  switching = true;
  /*Make a loop that will continue until
  no switching has been done:*/
  while (switching) {
    //start by saying: no switching is done:
    switching = false;
    rows = table.getElementsByTagName("TR");
    /*Loop through all table rows (except the
    first, which contains table headers):*/
    for (i = 1; i < (rows.length - 1); i++) {
      //start by saying there should be no switching:
      shouldSwitch = false;
      /*Get the two elements you want to compare,
      one from current row and one from the next:*/
      x = rows[i].getElementsByTagName("TD")[3];
      y = rows[i + 1].getElementsByTagName("TD")[3];
      //check if the two rows should switch place:
      if ((order == 'asc') && (Number(x.innerHTML.toLowerCase()) < Number(y.innerHTML.toLowerCase()))) {
        //if so, mark as a switch and break the loop:
        shouldSwitch = true;
        break;
      } else if ((order == 'desc') && (Number(x.innerHTML.toLowerCase()) > Number(y.innerHTML.toLowerCase()))) {
        //if so, mark as a switch and break the loop:
        shouldSwitch = true;
        break;
      }
    }
    if (shouldSwitch) {
      /*If a switch has been marked, make the switch
      and mark that a switch has been done:*/
      rows[i].parentNode.insertBefore(rows[i + 1], rows[i]);
      switching = true;
    }
  }
}
// --- SORT THE TABLE IN ASCENDING/DESCENDING ORDER, CREDIT: W3C. END ---

// --- SET TABLE ELEMENT ID'S START ---
function Set_Table_ID(probe) {
  // Change input id to match respective gene:
  var tochange = document.getElementById('myInput');
  tochange.id = probe + "_input";

  // Change table id to match respective gene:
  var tochangetwo = document.getElementById('myTable');
  tochangetwo.id = probe + "_table";

  // Change button sort id to match respective gene:
  var tochangethree = document.getElementById('myButton_asc');
  tochangethree.id = probe + "_button_asc";
  $("#" + tochangethree.id).change(() => {
    sortTable(probe, 'asc');
  });

  // Change button sort id to match respective gene:
  var tochangefour = document.getElementById('myButton_desc');
  tochangefour.id = probe + "_button_desc";
  $("#" + tochangefour.id).change(() => {
    sortTable(probe, 'desc');
  });
  // --- SET TABLE ELEMENT ID'S END ---
}


function Populate_Microarray_Table(affy_expression_dataset, affy_data, probe) {
  var labels = [];
  var labels = Object.keys(affy_data[0]);
  var affy_value =  Object.values(affy_data[0]);
  var description_table = document.getElementById("annotation_for_" + probe);
  
  for(let i=4; i<labels.length; i++) {    // Retrieve affy values for graph data from susbset of array.
    var temp = [];
    temp.push(labels[i], affy_value[i]);
    affy_expression_dataset.push(temp);
  }

  for(let i = 0; i<3; i++){ // Fill description table with sql data from subset of array.
    description_table.rows[i].children[1].innerHTML = affy_value[i];
    description_table.rows[i].children[0].style.background = "#4286f4";
    description_table.rows[i].children[0].style.color = "white";
  }


  let desc_array = affy_value[3].split('/');  // Replace description / with , for better presentation.
  desc_array.shift();
  description_table.rows[3].children[1].innerHTML = desc_array[6];
  description_table.rows[3].children[0].style.background = "#4286f4";
  description_table.rows[3].children[0].style.color = "white";

  return affy_expression_dataset;
}
