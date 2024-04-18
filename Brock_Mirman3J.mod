// Declare endogenous variables
var c k z y invest;

// Declare exogenous shock
varexo e_z;

// Declare parameters
parameters beta delta alpha rho sigma;

// Assign parameter values
beta = 0.95;   // Discount factor
delta = .05;  // Depreciation rate
alpha = 0.36;  // Capital share in production
rho = 0.9;     // Autoregressive coefficient for the shock
sigma = 3;   // CRRA coefficient (adjust as needed)

// Model equations
model;
  // Euler equation
  c^(-sigma) = beta * c(+1)^(-sigma) * (alpha *exp(z) * k^(alpha - 1) + 1 - delta);

  // Capital accumulation equation
  k = (1 - delta) * k(-1) + exp(z) * k(-1)^alpha - c;

  // Stochastic process for the shock
  z = rho * z(-1) + e_z;
  y = exp(z) * k(-1)^alpha;
  invest = y - c;
end;

//options_.solve_toler = 1e-8; 
//options_.solve_maxit = 10000;
//steady;
// Initial values (steady state)
initval;
  // Calculate steady state values
  k = ((1 - beta * (1 - delta)) / (beta * alpha)) ^ (1 / (alpha - 1));
  c = k^alpha - delta * k;
  z = 0;
  y = y^alpha;
  invest = y - c;
end;
steady;
check;

// Stochastic setup
shocks;
var e_z; stderr 0.05;
end;

// Simulate
stoch_simul(order=1, periods=1000, irf=0);
CSV.write("BrockMirman.csv", context.results.model_results[1].simulations[1].data)

//varobs y;
//calib_smoother(datafile='BrockMirman.csv', diffuse_filter, filtered_vars);

