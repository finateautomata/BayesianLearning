## @knitr sim_adoption
library(tidyverse)

gen_y_discrete <- function(n,lambda){
  y <- as.integer(runif(n) < lambda)
}

Nsim <- 10000
lambda <- runif(Nsim,min=0.01,max=0.50)
n <- 100

simDF <-  NULL

for (i in 1:Nsim){
  y <- gen_y_discrete(n,lambda[i])
  simDF <- bind_rows(simDF,
                     data.frame(y = sum(y), lambda = lambda[i]))
}

simDF %>%
  ggplot(aes(x = y, y = lambda)) + geom_jitter(alpha=0.3) + 
  labs(x = "Number of Adopters in Sample",
       y = "True Adoptation Rate")


## @knitr hist_cond_data

# One survey sample says 20 adopters - what are plausible values of p?
simDF %>%
  filter(y == 20) %>%
  ggplot(aes(x = lambda,y=..density..)) + geom_histogram() + 
  geom_density(alpha=.2, fill="#FF6666") + xlim(0,0.40) + 
  labs(title = 'Distribution of lambda conditional on y=20 adopters')

# One survey sample says 10 adopters - what are plausible values of lambda?
simDF %>%
  filter(y == 10) %>%
  ggplot(aes(x = lambda,y=..density..)) + geom_histogram() + 
  geom_density(alpha=.2, fill="#FF6666") + xlim(0,0.40) + 
  labs(title = 'Distribution of lambda conditional on y=10 adopters')

## @knitr calcProb

yObs <- 20
pUnder <- 0.3

tmp <- simDF %>%
  filter(y == yObs)

probP <- sum(tmp$lambda > pUnder)/nrow(tmp)

## @knitr hist_estimator

# True lambda=0.20

lambda0 = 0.20 
lambdahat <- rep(0,Nsim)

for (i in 1:Nsim){
  y <- gen_y_discrete(n,lambda0)
  lambdahat[i] <- sum(y)/n
}

data.frame(lambdahat = lambdahat) %>%
  ggplot(aes(x = lambdahat,y=..density..)) + geom_histogram() + 
  geom_density(alpha=.2, fill="#FF6666") + xlim(0,0.40) + 
  labs(title = 'Distribution of Estimator when lambda=0.2')

## @knitr plotBetaPost

N1 <- 20
N <- 100
lambdaSeq <- seq(0.01,0.99,0.01)

data.frame(lambda = lambdaSeq, posterior = dbeta(lambdaSeq,N1+1,N-N1+1), prior = 1) %>%
  pivot_longer(names_to = "Distribution",values_to = "Value", -lambda) %>%
  ggplot(aes(x = lambda,y = Value, color = Distribution)) + geom_line() + xlim(0,0.5) +
  labs(title = 'Posterior Distribution of lambda, N1 = 20')

## @knitr bayesAB

# data
yA <- 317
nA <- 10000
yB <- 152
nB <- 5000

# prior 
a0 <- 1
b0 <- 20

# posterior draws of CTR for A and B 
nSim <- 10000
lambdaAPost <- rbeta(nSim,yA + a0,nA - yA + b0)
lambdaBPost <- rbeta(nSim,yB + a0,nB - yB + b0)
deltaPost <- lambdaAPost - lambdaBPost

ProbDeltaPos <- sum(deltaPost > 0)/nSim


postDF <- data.frame(lambdaA = lambdaAPost,
                     lambdaB = lambdaBPost,
                     delta = deltaPost
                     ) %>%
  rowid_to_column(var = 'draw')

postDF %>%
  pivot_longer(names_to = 'Parameter',values_to = 'Value', -draw) %>%
  ggplot(aes(x=Value,fill=Parameter)) + geom_density(alpha=0.3) + 
  geom_vline(aes(xintercept=0)) + 
  labs(title = 'Posterior Distributions')


## @knitr posteriorLoss

lossF <- function(choice,lA,lB){
  
  if (choice=='A'){
    loss <- mean((lB - lA)*(lB > lA))
  } else {
    loss <- mean((lA - lB)*(lA > lB))
  }
  
  return(loss)
}

lossA <- lossF('A',lambdaAPost,lambdaBPost)
lossB <- lossF('B',lambdaAPost,lambdaBPost)



## @knitr plotMuPost

#mu <- 1.0
sigma <- 1.0
#N <- 100
# y <- rnorm(N,mu,sigma)
# saveRDS(y, 'data/y.rds')

y <- read_rds('data/y.rds')

postDens <- function(y,sigma,mu0,sigma0){
  
  N <- length(y)
  muMin <- 0.0
  muMax <- 2.0
  muSeq <- seq(muMin,muMax,0.01)
  
  precPost <- (N/(sigma^2)) + (1.0/(sigma0^2))
  mnPost <- ((sum(y)/(sigma^2)) + (mu0/(sigma0^2)))/precPost
  
  return( data.frame(mu = muSeq, 
                     density = dnorm(muSeq,mean = mnPost,sd = sqrt(1.0/precPost)))
  )
  
}

df <- bind_rows(
  postDens(y,sigma,0.0,10.0) %>%
    mutate(mu0 = 0.0,sigma0 = 10.0),
  postDens(y,sigma,0.0,0.5) %>%
    mutate(mu0 = 0.0,sigma0 = 0.5),
  postDens(y,sigma,0.0,0.25) %>%
    mutate(mu0 = 0.0,sigma0 = 0.25),
  postDens(y,sigma,0.0,0.1) %>%
    mutate(mu0 = 0.0,sigma0 = 0.1)
)



df %>%
  mutate(sigma0 = factor(sigma0)) %>%
  ggplot(aes(x = mu, y = density, color = sigma0)) + geom_line() + geom_vline(aes(xintercept = 1.0)) +
  labs(title = 'Posterior for mu', subtitle = 'N = 100, sigma = 1.0, mu0 = 0.0')







