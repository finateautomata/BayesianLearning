data {
  int<lower=0> nObs;                         // number of rows in full data 
  int<lower=0> nStores;                      // number of stores 
  int<lower=1,upper=nStores> storeID[nObs];  // store index for each row
  int<lower=1> dimZ;                         // dimension of store level predictors
  vector[nObs] lp;                           // log price
  vector[nObs] ls;                           // log sales volume
  row_vector[dimZ] Z[nStores];               // group predictors
}

parameters {
  corr_matrix[2] Omega;        // prior correlation
  vector<lower=0>[2] tau;      // prior scale
  vector[2] theta[nStores];    // indiv coeffs by group
  real<lower=0> sigma;         // prediction error scale
  matrix[dimZ, 2] gamma;       // group coeffs
}

model {
  tau ~ cauchy(0, 2.5);
  Omega ~ lkj_corr(2);
  to_vector(gamma) ~ normal(0, 5);
  
  {
    row_vector[2] z_gamma[nStores];
    for (s in 1:nStores) z_gamma[s] = Z[s] * gamma;
    theta ~ multi_normal(z_gamma, quad_form_diag(Omega, tau));
  }

  sigma ~ cauchy(0, 2.5);
  
  for (n in 1:nObs)
    ls[n] ~ normal(theta[storeID[n]][1] + lp[n]*theta[storeID[n]][2], sigma);
}
