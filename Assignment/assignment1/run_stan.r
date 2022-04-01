## run stan model 
library(tidyverse)
library(rstan)

df <- read_csv('data.csv')

bilModel <- stan_model(file='MLM_binomial.stan')

Niter <-  1500
Nwarmup <- 1000
Nchains <- 4

fit <- sampling(bilModel,
                data=list(nObs = nrow(df),
                          nGroups = max(df$group),
                          groupID = df$group,
                          y = df$y),
                iter = Niter,
                warmup = Nwarmup,
                chains = Nchains,
                verbose = TRUE)
## save draws

thePars <- rstan::extract(fit,pars = c('theta','mu','sigma'))

thetaDraws <- as.data.frame(thePars$theta) %>%
  rowid_to_column(var = 'draw') %>%
  pivot_longer(-draw,names_to = 'group',values_to = 'value') %>%
  separate(group,into = c('tmp','group'),sep=1) %>%
  select(-tmp) %>%
  mutate(group = as.integer(group))

muSigDraws <- data.frame(mu = thePars$mu, sigma = thePars$sigma, draw = 1:length(thePars$mu))

write_csv(thetaDraws, path = 'thetaDrawsM1.csv')  
write_csv(muSigDraws, path = 'muSigmaDrawsM1.csv')  


