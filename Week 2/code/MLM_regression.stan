data {
  int<lower=0> nObs;                         // number of rows in full data 
  int<lower=0> nStores;                      // number of stores 
  int<lower=1,upper=nStores> storeID[nObs];  // store index for each row
  vector[nObs] lp;                           // log price
  vector[nObs] ls;                              // log sales volume
}

parameters {
  corr_matrix[2] Omega;        // prior correlation
  vector<lower=0>[2] tau;      // prior scale
  matrix[nStores,2] theta;     // store specific coefficients 
  real<lower=0> sigma;         // prediction error scale
  row_vector[2] mu;            // mean coefficients
}

model {
  tau ~ cauchy(0, 2.5);
  Omega ~ lkj_corr(2);
  mu ~ normal(0,5);
  
  for (s in 1:nStores) {
    theta[s,] ~ multi_normal(mu, quad_form_diag(Omega, tau));
  }

  sigma ~ cauchy(0, 2.5);
  
  for (n in 1:nObs)
    ls[n] ~ normal(theta[storeID[n]][1] + lp[n]*theta[storeID[n]][2], sigma);
}
