data {
  int<lower=0> nObs;                         // number of rows in full data 
  int<lower=0> nStores;                      // number of stores 
  int<lower=1,upper=nStores> storeID[nObs];  // store index for each row
  vector[nObs] lp;                           // log price
  vector[nObs] ls;                              // log sales volume
}

parameters {
  matrix[nStores,2] theta;     // store specific coefficients 
  real<lower=0> sigma;         // prediction error scale
}

model {
  
  theta[,1] ~ normal(0, 10);
  theta[,2] ~ normal(0, 10);
  sigma ~ cauchy(0, 2.5);
  
  for (n in 1:nObs)
    ls[n] ~ normal(theta[storeID[n]][1] + lp[n]*theta[storeID[n]][2], sigma);
}
