data {
  int<lower=0> nObs;                         // number of rows in full data 
  int<lower=0> nUsers;                       // number of users
  int<lower=1,upper=nUsers> userID[nObs];    // user index for each row
  vector[nObs] y;                            // log amount
}

parameters {
  real<lower=0> sigma;         // sd alpha
  real mu;                     // mean alpha
  vector[nUsers] alpha;        // user effects
  real<lower=0> sigma_y;       // sd data
}

model {
  sigma ~ cauchy(0, 2.5);
  mu ~ normal(0,5);
  alpha ~ normal(mu, sigma);
  sigma_y ~ cauchy(0, 2.5);

  y ~ normal(alpha[userID], sigma_y);
}
