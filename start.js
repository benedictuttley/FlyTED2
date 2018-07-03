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

app.get('/', (req, res) => {res.sendFile('/home/benedict/Desktop/index.html')});


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

//
// res.render('hello', {
//     query_result: [
//       "Yehuda Katz",
//       "Alan Johnson",
//       "Charles Jolley"
//     ]});
//   });

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

  console.log("General Query received");

  fetchAll((err, result) => {
    if (err) {
      throw err
    } else {
      res.send(result)
    }
  });
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
  }, req.body.xAxis);
});

// List for probe query post from HTML form:
app.post('/y_query', (req, res) => {

  console.log("Y-coordinate Query received");

  fetchX((err, result) => {
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


// [1] Fetch all entries in table and list them:
function fetchAll(callback){
  conn.query("SELECT * FROM Demo", (err, result) => {
    if (err)
     throw err;
    else {
      callback(null, result)
    }
  });
}

// [2] Fetch all entries with a certain Probe name and list them:
function fetchProbe(callback, probeName){
  console.log(probeName);
  conn.query("SELECT * FROM Demo WHERE Probe = " + "'" + probeName + "'" , (err, result) => {
    if (err)
     throw err;
    else {
      callback(null, result)
    }
  });
}

// [3] Fetch all entries that have a particular x coordinate greater than a given number:
function fetchX(callback, xAxis){
  console.log(xAxis);
  conn.query("SELECT * FROM Demo WHERE Xcoordinate >= " + "'" + xAxis + "'", (err, result) => {
    if (err)
     throw err;
    else {
      callback(null, result)
    }
  });
}

// [4] Fetch all entries that have a particular y coordinate greater than a given number:
function fetchX(callback, yAxis){
  console.log(yAxis);
  conn.query("SELECT * FROM Demo WHERE Ycoordinate >= " + "'" + yAxis + "'", (err, result) => {
    if (err)
     throw err;
    else {
      callback(null, result)
    }
  });
}

// [4] Fetch all entries that have a particular y coordinate greater than a given number:
function fetchGenA(callback, genA){
  console.log(genA);
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
