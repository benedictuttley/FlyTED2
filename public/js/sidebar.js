// --- SIDEBAR TOGGLE SCRIPT START ---
function toggleSide() {
  let side = document.getElementById("mySidenav");
  if (side.style.width === "250px") { side.style.width = "0px"} // Close the sidenav if it was previously open.
  else {
    side.style.width = "250px" // Open the sidenav if it was previously closed.
    side.style.borderRight = "2px solid white";
  }
}

function closeNav() {
  document.getElementById("mySidenav").style.width = "0";
}

// --- SIDEBAR TOGGLE SCRIPT END ---
