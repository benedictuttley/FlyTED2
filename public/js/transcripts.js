// --- EDIT STYLE OF SEQUENCES TO SHOW 5' AND 3' PROBES START ---
function show_probes(Transcript_ID, Three_Prime_Sequence, Five_Prime_Sequence) {

  let test = document.getElementById(Transcript_ID);
  let three = test.innerHTML.replace(Three_Prime_Sequence, "");
  let res = three.replace(Five_Prime_Sequence, "").trim();

  let final =
    '<span style="color: #660066; text-decoration: underline overline; font-weight: bold">' +
    Five_Prime_Sequence +
    '</span>' + '<span>' + res + '</span>' +
    '<span style="color: #009933; text-decoration: underline overline; font-weight: bold;">' +
    Three_Prime_Sequence + '</span>';

  document.getElementById(Transcript_ID).innerHTML = final;
  document.getElementById(Transcript_ID + "_container").style.display = "none"; // Initially have transcripts hiddden.
}

// --- EDIT STYLE OF SEQUENCES TO SHOW 5' AND 3' PROBES END ---

// --- TOGGLE DISPOLAY OF TRANSCRIPT CONTENTS START ---

function toggle_transcript_view(transcript_id) {
  let transcipt_view = document.getElementById(transcript_id);
  if (transcipt_view.style.display === "none") {
    transcipt_view.style.display = "block";
  } else {
    transcipt_view.style.display = "none";
  }
}
// --- TOGGLE DISPOLAY OF TRANSCRIPT CONTENTS END ---


function Create_Transcript_Header() {
  let probe_header = document.getElementsByTagName("div");
  let transcript_header = document.createElement("h2");
  let text = document.createTextNode("TRANSCRIPTS");
  transcript_header.style.color = "#d4374a";
  transcript_header.appendChild(text);
  probe_header[probe_header.length - 1].appendChild(transcript_header); // APPEND TRANSCRIPTS HEADER TO RESPECTIVE GENE HEADER.
}

function Create_Transcript_Tab(Transcript_ID) {
  let transcript_button = document.getElementById('temp');
  let id = Transcript_ID + "_tab"; // Set the transcript tab equal to the transcript ID.
  transcript_button.setAttribute("id", id);
  transcript_button.classList.add('tab');
  transcript_button.innerHTML =
    ">> VIEW TRANSCRIPT VIEW INFORMATION FOR GENE: " +
    Transcript_ID;

  transcript_button.onclick = function() { // When tab is clicked toggle view of the respective transcipt information.
    toggle_transcript_view(Transcript_ID + "_container");
  };
}
