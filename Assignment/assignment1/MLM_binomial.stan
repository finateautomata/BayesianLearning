data {
  int<lower=0> nObs;                         // number of rows in full data 
  int<lower=0> nGroups;                      // number of groups
  int<lower=1,upper=nGroups> groupID[nObs];  // group index for each row
  int y[nObs];                            // outcomes 
}

parameters {
  vector[nGroups] theta;     // group specific coefficients 
  real<lower=0> sigma;       // prediction error scale
  real  mu;                  // mean coefficients
}

model {
  sigma ~ cauchy(0, 2.5);
  mu ~ normal(0,5);
  theta ~ normal(mu,sigma);
  
  y ~ bernoulli_logit(theta[groupID]);
}
