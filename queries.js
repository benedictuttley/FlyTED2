exports.createQuery = function(probeName) {

  var expression_tissue_query = { // Tissue expresssion query
    "name": "Gene_FlyAtlas",
    "title": "Gene --> FlyAtlas expression data",
    "comment": "06.06.07 updated to work from gene class - Philip",
    "description": "For a given D. melanogaster gene, show expression data from FlyAtlas.",
    "constraintLogic": "A and B and C and D and E and F",
    "from": "Gene",
    "select": [
      "secondaryIdentifier",
      "symbol",
      "microArrayResults.mRNASignal",
      "microArrayResults.mRNASignalSEM",
      "microArrayResults.presentCall",
      "microArrayResults.enrichment",
      "microArrayResults.affyCall",
      "microArrayResults.dataSets.name",
      "microArrayResults.tissue.name",
      "primaryIdentifier"
    ],
    "orderBy": [{
      "path": "microArrayResults.tissue.name",
      "direction": "ASC"
    }],
    "where": [{
        "path": "organism.name",
        "op": "=",
        "value": "Drosophila melanogaster",
        "code": "B",
        "editable": false,
        "switched": "LOCKED",
        "switchable": false
      },
      {
        "path": "microArrayResults",
        "type": "FlyAtlasResult"
      },
      {
        "path": "Gene",
        "op": "LOOKUP",
        "value": probeName,
        "code": "A",
        "editable": true,
        "switched": "LOCKED",
        "switchable": false
      },
      {
        "path": "microArrayResults.affyCall",
        "op": "NONE OF",
        "values": [
          "None"
        ],
        "code": "C",
        "editable": true,
        "switched": "ON",
        "switchable": true
      }
    ]
  }
  return expression_tissue_query;
}
