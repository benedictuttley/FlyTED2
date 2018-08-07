// --- EDIT STYLE OF SEQUENCES TO SHOW 5' AND 3' PROBES START ---
function show_probes(Transcript_ID, Three_Prime_Sequence, Five_Prime_Sequence) {
  let test = document.getElementById(Transcript_ID);
  let three = test.innerHTML.replace(Three_Prime_Sequence, "");
  let res = three.replace(Five_Prime_Sequence, "").trim();

  let final =
    '<span style="color: #660066; text-decoration: underline overline; font-weight: bold">' +
    "{{[5_Prime_Sequence]}}" +
    '</span>' + '<span>' + res + '</span>' +
    '<span style="color: #009933; text-decoration: underline overline; font-weight: bold;">' +
    "{{[3_Prime_Sequence]}}" + '</span>';

  document.getElementById(Transcript_ID).innerHTML = final;
  document.getElementById(Transcript_ID + "_container").style.display = "none"; // Initially have transcripts hiddden.
}

// --- EDIT STYLE OF SEQUENCES TO SHOW 5' AND 3' PROBES END ---
