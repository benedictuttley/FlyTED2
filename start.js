// Note: Currently the primary key consists of: Slide Name, Genotype A, Genotype B, X-coordinate and
// Y-coordinate.
const express = require('express')
const app = express();
const mysql = require('mysql');
const bodyParser = require("body-parser");

// Include handlbars
const handlebars = require('express-handlebars');
const Handlebars = require('handlebars');
app.engine('handlebars', handlebars({defaultLayout: 'main'}));
app.set('view engine', 'handlebars');

app.use(bodyParser.urlencoded({
    extended: true
}));

app.use(bodyParser.json());

app.post("/test", function (req, res) {
    console.log(req.body.name)
    res.send('Understood');
});

app.get('/', (req, res) => {res.sendFile('/home/benedict/Desktop/FlyTED2/index.html')});


// Demo handlebar helper function:
Handlebars.registerHelper('with', (row) => {
  return JSON.stringify(row);
});


app.get('/htest', (req,res) => {

  fetchAll((err, result) => {
    if (err) {
      throw err
    } else {
      res.render('hello', {query_result: result});
    }
  });
});

app.listen(3000, () => console.log('App active on port:3000'))

// Create mysql connection to existing FlyTed database:
var conn = mysql.createConnection({
    host: 'localhost',
    user: 'pmauser',
    password: '[*3G4vlSp7',
    database: 'FlyTed2'
});

// Test connection:
conn.connect(err => {
    if (err) throw err;
    console.log("Mysql Connection Established");
});

// Query Database (simple) - Actual Function:
function test(callback) {
    conn.query('SELECT * FROM table_name_for_df ', (err, result) => {
        if (err) throw err;
        if (err)
            callback(err, null);
        else
            callback(null, result);
    });
}

// Attempt to send FlyTed2 image to client and render on screen:
// Need to:
// Fetch from DB
// Test Vue with db


// Listen for general query post from HTML form:
app.post('/query', (req, res) => {

  fetchAll((err, result) => {
    if (err) {
      throw err
    } else {
      res.send(result)
    }
  }, req.body.Group);
});


// List for probe query post from HTML form:
app.post('/probe_query', (req, res) => {

  console.log("Probe Query received");

  fetchProbe((err, result) => {
    if (err) {
      throw err
    } else {
      res.send(result)
    }
  }, req.body.Probe);
});

// List for probe query post from HTML form:
app.post('/x_query', (req, res) => {

  console.log("X-coordinate Query received");

  fetchX((err, result) => {
    if (err) {
      throw err
    } else {
      res.send(result)
    }
  }, req.body.x_min, req.body.x_max);
});

// List for probe query post from HTML form:
app.post('/y_query', (req, res) => {

  console.log("Y-coordinate Query received");

  fetchY((err, result) => {
    if (err) {
      throw err
    } else {
      res.send(result)
    }
  }, req.body.yAxis);
});

// List for genotype b query post from HTML form:
app.post('/gena_query', (req, res) => {

  console.log("Genotype A Query received");

  fetchGenA((err, result) => {
    if (err) {
      throw err
    } else {
      res.send(result)
    }
  }, req.body.genA);
});


// List for genotype b query post from HTML form:
app.post('/genb_query', (req, res) => {

  console.log("Genotype B Query received");

  fetchGenB((err, result) => {
    if (err) {
      throw err
    } else {
      res.send(result)
    }
  }, req.body.genB);
});


// TODO: [1] Fetch all entries in table and list them:
function fetchAll(callback, Group){
  var group_by;
  Group.length == "None" ? group_by = (" GROUP BY " + Group) : group_by = "";

  conn.query("SELECT * FROM Demo" + group_by, (err, result) => {
    if (err){
      throw err;
    }
    else {
      callback(null, result);
    }
  });
}

// TODO: [1] Fetch all entries and group them by probe name:
function group_by_probe(){
  conn.query("SELECT * FROM Demo GROUP BY Probe", (err, result) => {
    if(err){
      throw err;
    }
    else {
      callback(null,result);
    }
  });
}

// TODO: Fetch all entries and group them by Genotype A:
function group_by_genA() {
  conn.query("SELECT * FROM Demo GROUP BY GenotypeA", (err, result) => {
    if(err){
      throw err;
    }
    else {
      callback(null, result);
    }
  });
}

// TODO: Fetch all entries and group them by Genotype B:
function group_by_genB() {
  conn.query("SELECT * FROM Demo GROUP BY GENOTYPE B", (err, result) => {
    if(err){
      throw err;
    } else {
      callback(null, result);
    }
  });
}


// [2] Fetch all entries with a certain Probe name and list them:
// [2.a] Select enries corresponding to multiple probes:
function fetchProbe(callback, probeName){
  var formatted_probes = ""
  probeName = probeName.split(", ");
  probeName.forEach(probe => formatted_probes += ("'" + probe + "' "));
  formatted_probes = news.trim().replace(/\s+/g, ', ');

  conn.query("SELECT * FROM Demo WHERE Probe IN (" + formatted_probes + ")" , (err, result) => {
    if (err)
     throw err;
    else {
      callback(null, result)
    }
  });
}

// [3] Fetch all entries that have a particular x-coordinate greater than a given number:
function fetchX(callback, xMin, xMax){

  // Check if min x-value or max x-value input is empty, if true then permit any lower/upper bound:
  xMin.length == 0 ? xMin = 0 : xMin = xMin;
  xMax.length == 0 ? xMax = "(SELECT MAX(Xcoordinate) FROM Demo)": xMax = xMax;

  conn.query("SELECT * FROM Demo WHERE Xcoordinate >= " + "'" + xMin + "'" + " AND Xcoordinate <= " + xMax , (err, result) => {
    if (err)
      throw err;
    else {
      callback(null, result)
    }
  });
}


// [4] Fetch all entries that have a particular y-coordinate greater than a given number:
function fetchY(callback, yMin, yMax){

  // Check if min y-value or max y-value input is empty, if true then permit any lower/upper bound:
  yMin.length == 0 ? yMin = 0 : yMin = y_Min;
  yMax.length == 0 ? yMax = "SELECT MAX(Ycoordinate) FROM Demo" : yMax = y_Max;

  conn.query("SELECT * FROM Demo WHERE Ycoordinate >= " + "'" + yMin + "'" + " AND Ycoordinate <= " +  yMax , (err, result) => {
    if (err)
     throw err;
    else {
      callback(null, result)
    }
  });
}


// [4] Fetch all entries that have a particular y coordinate greater than a given number:
function fetchGenA(callback, genA){
  co
  conn.query("SELECT * FROM Demo WHERE GenotypeA = " + "'" + genA + "'", (err, result) => {
    if (err)
     throw err;
    else {
      callback(null, result)
    }
  });
}

// [4] Fetch all entries that have a particular y coordinate greater than a given number:
function fetchGenB(callback, genB){
  console.log(genB);
  conn.query("SELECT * FROM Demo WHERE GenotypeB = " + "'" + genB + "'", (err, result) => {
    if (err)
     throw err;
    else {
      callback(null, result)
    }
  });
}


// Compound query builder:
// Variables to consider: probe, slide, Genotype A and B, Date, Coordinates.

// // TODO: Coordinates Options
// --------------------

// Options for coordinates: The range, from n to p, n to max, as a slider:
// If unchecked then coordinate can be any value
