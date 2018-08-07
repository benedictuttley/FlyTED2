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
