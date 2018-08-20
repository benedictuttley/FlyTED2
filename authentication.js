const passport = require('passport');
const LocalStrategy = require('passport-local').Strategy;
const expressSession = require('express-session');
const flash = require('express-flash-messages');

// IN MEMORY USER
const Users = {
  Helen: {
    username: 'Helen',
    password: 'Cardiff_2018'
  }
}
// HASH THE IM MEMORY USER'S PASSWORD
bcrypt.hash(Users.Helen.password, 10, (err, hash) => {
  Users.Helen.password = hash;
});

// When user attempts login, authenticate with Passport.js
passport.use("login", new LocalStrategy({
    passReqToCallback: true // Allow Passport.js to use the request object
  },
  (req, username, password, done) => {

    findUser(username, (err, user) => {
      if (err)
        return done(err);
      if (!user) {
        return done(null, false)
      }

      bcrypt.compare(password, user.password, (err, res) => {
        if(res){
          console.log("THE PASSWORDS MATCH!");
          return done(null, user); // On correct credentials
          }
        else{
          console.log("THE PASSWORDS DONT MATCH");
          return done(null, false, req.flash('error', 'This is a flash message using the express-flash module.'));
        }
      });
    });
  }));

passport.serializeUser((user, done) => {  // Create the user for the session.
  var liteUser = createLiteUser(user);
  done(null, liteUser);
});

passport.deserializeUser((user, done) => { // Remove the user from the session.
  done(null, user);
});

function createLiteUser(user){
  let liteUser = {
    email: "test"
  }
  return liteUser;
}

function findUser(username, done) { // Find the user object by username.
  done(null, Users[username]);
}

// Check if user is authenticated or not to determine which admin page view to serve.
function ensureAuthenticated(req, res, next) {
  if (req.isAuthenticated()) {
    return next();
  } else{
    // Not authenticated then re-serve the admin page aswell as including a flash description.
    req.flash('error', 'You must be logged in to view this resource.');
    var flashMessages = res.locals.getMessages(); // Fetch errors sored in flash memory.
    res.render('admin', {
      errors: flashMessages
    });
  }
}


module.exports = passport;s
