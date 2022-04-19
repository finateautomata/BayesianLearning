## @knitr loadLibraries
library(tidyverse)


## @knitr exampleMC

rho <- 0.75
sigma <- 1.0 

nSim <- 5000
eps <- rnorm(nSim, mean = 0, sd = sigma)

# initialize array and set beginning
df <- data.frame(x = array(0,nSim))
df$step <- 1

df$x[1] <- 2.5

for (t in 2:nSim){
  df$x[t] = rho*df$x[t-1] + eps[t]
  df$step[t] <- t
}

init.pl <- df %>%
  slice(1:20) %>%
  ggplot(aes(x = step, y = x)) + geom_line() + 
  labs(title = 'First 20 Steps')

all.pl <- df %>%
  ggplot(aes(x = step, y = x)) + geom_line() + 
  labs(title = paste0('All ',nSim,' Steps'))

StationaryMean <- 0
StationarySD <- sqrt(sigma^2/(1-rho^2))

nStart <- 1001
nSteps <- 4000
simMean <- mean(df$x[nStart:(nStart + nSteps - 1)])
simSD <- sd(df$x[nStart:(nStart + nSteps - 1)])


## @knitr exampleMetropolis 
library(mvtnorm)

# target distribution 
target <- function(x){

    mu1 <- c(-0.5,0)
    sigma1 <- array(c(.75,0.25,0.25,.75),c(2,2))
    
    mu2 <- c(0.25,1.5)
    sigma2 <- array(c(0.5,-0.25,-0.25,0.5),c(2,2))

    l <- 0.6
        
    d <- l*dmvnorm(x, mean = mu1, sigma = sigma1) + 
      (1.0-l)*dmvnorm(x, mean = mu2, sigma = sigma2)

}

# plot target 
points <- data.frame(expand.grid(x1 = seq(-3,3,.01),x2 = seq(-3,3,.01))) 
points$d <- target(points[,1:2])

# ggplot contour
ggpl <- points %>%
  ggplot() +
  geom_tile(aes(x=x1,y=x2,fill=d)) +
  geom_contour(aes(x=x1,y=x2,z=d),color="black") +  
  scale_fill_gradientn("Z",colours = terrain.colors(10)) 

## 
## rayshader version 
##
## note: takes a bit of time to render 

#library(rayshader)
#plot_gg(ggpl, multicore = FALSE, raytrace = TRUE, width = 7, height = 4,
#        scale = 300, windowsize = c(1400, 866), zoom = 0.6, phi = 30, theta = 30)
#render_snapshot(clear = TRUE)

## @knitr sampleExampleMetropolis 

# jump factor (scale to get good acceptance rates)
alpha <- 2.0

nSim <- 25000
chain <- array(0,c(nSim,2))
accept <- array(0,c(nSim))

for (i in 2:nSim){
  
  ## propose new move 
  e <- rnorm(2, sd = alpha)
  x_star <- chain[i-1,] + e
  
  ## calculate accept rate 
  dCurrent <- target(chain[i-1,])
  dNew <- target(x_star)
  
  r <- dNew/dCurrent

  if (runif(1) < r){
    chain[i,] <- x_star
    accept[i] <- 1
  } else {
    chain[i,] <- chain[i-1,]
  }
  
}

chainDF <- data.frame(chain)
burnIn <- 5000
saveRDS(chainDF, file = 'data/exampleMetropolisChain.rds')
saveRDS(accept, file = 'data/accept.rds')


## @knitr sampleExampleMetropolis2

chainDF <- read_rds('data/exampleMetropolisChain.rds')
accept <- read_rds('data/accept.rds')

nSim <- nrow(chainDF)
burnIn <- 5000
acceptRate <- mean(accept[(burnIn+1):nSim])

ggpl2 <- ggpl  + 
  geom_point(data=chainDF[(burnIn+1):nSim,],aes(x = X1,y=X2),alpha=0.3, size=0.3) + 
  labs(title = 'Last 20,000 Iterations')

## @knitr metropolisLogitModelData 

# generate synthetic data for logit classifier
N <- 1000
b <- c(0.0,0.5,1.0,-0.5,0.0)
X <-  cbind(array(1,c(N,1)), array(rnorm(N*4),c(N,4)) )
Xb <- X %*% b
pr <- exp(Xb)/(1.0 + exp(Xb))
y <- runif(N) < pr 

df <- data.frame(y,X)
saveRDS(df, file = 'data/logitData.rds')
write_csv(df, path = 'data/logitData.csv')

## @knitr metropolisLogitModel

# function returning log-posterior for a given beta
logPost <- function(y,X,tau,beta){
  
  Xb <- X %*% beta
  
  # log-likelihood
  ll <- sum(Xb[y]) - sum(log(1.0 + exp(Xb)))
  
  # log-posterior 
  lp <- ll - 0.5*tau*sum(beta^2)
  
  return(lp)
}

# read data 
df <- read_rds('data/logitData.rds')
y <- df$y
X <- as.matrix(df[,2:6])
K <- ncol(X)

# set prior precision = 1.0/(prior variance)
tau <- 0.2 

# set jump factor 
alpha <- 0.08

nSim <- 15000
chain <- array(0,c(nSim,K))
accept <- array(0,c(nSim))


for (i in 2:nSim){
  
  ## propose new move 
  x_star <- chain[i-1,] + rnorm(K, sd = alpha)
  
  ## calculate accept rate 
  lpCurrent <- logPost(y,X,tau,chain[i-1,])
  lpNew <- logPost(y,X,tau,x_star)
  logR <- lpNew - lpCurrent
  
  ## accept/reject step
  if (log(runif(1)) < logR){
    chain[i,] <- x_star
    accept[i] <- 1
  } else {
    chain[i,] <- chain[i-1,]
  }
  
}

chainDF <- data.frame(chain)
nSim <- nrow(chainDF)
burnIn <- 5000
acceptRate <- mean(accept[(burnIn+1):nSim])

saveRDS(chainDF, file = 'data/logitMetropolis.rds')
saveRDS(accept, file = 'data/logitAccept.rds')

## @knitr plotLogitMetropolis 

chainDF <- read_rds('data/logitMetropolis.rds')
accept <- read_rds('data/logitAccept.rds')
burnIn <- 5000

trueValues <- data.frame(value = c(0.0,0.5,1.0,-0.5,0.0))
trueValues$parameter <- paste0('X',c(1:5))

chainDF %>%
  rowid_to_column(var = 'draw') %>%
  filter(draw > burnIn) %>%
  pivot_longer(-draw,names_to = 'parameter', values_to='value') %>%
  ggplot(aes(x = draw,y= value, color = parameter)) + geom_line() + 
  geom_hline(data=trueValues,aes(yintercept = value)) + 
  facet_wrap(~parameter) + 
  theme(legend.position = "none") + 
  labs(x = 'Iteration',
       title = '10,000 Last Iterations',
       subtitle = 'True parameter value indicated by horizontal line')

## @knitr generateDataGibbs1

# first set true parameters that we are trying to learn about
mu <- 1.0     # mean
tau <- 1.0    # precision = 1/variance 

# set sample size and generate data
N <- 250
df <- data.frame(X = rnorm(N,mean = mu, sd = sqrt(1.0/tau)))
write_csv(df, path = 'data/dataGibbs1.csv')


## @knitr runGibbs1

# read data and pre-compute sum
df <- read_csv('data/dataGibbs1.csv')
N <- nrow(df)
sum_x <- sum(df$X)

# set prior
mu0 <- 0.0
tau0 <- 0.2
a0 <- 1.5
b0 <- 0.5

# start sampler 
nIter <- 3000
gibbsDF <- data.frame(mu = array(0,nIter), tau = array(0,nIter))
gibbsDF$tau[1] <- 0.5  # initialization 

for (i in 2:nIter){
  ## mu|tau
  mu_var <- 1.0/(gibbsDF$tau[i-1]*N + tau0)
  mu_mean <- mu_var*(gibbsDF$tau[i-1]*sum_x + tau0*mu0)
  gibbsDF$mu[i] <- rnorm(1,mean = mu_mean, sd = sqrt(mu_var))
  
  ## tau|mu
  ag <- a0 + 0.5*N
  bg <- b0 + 0.5*sum((df$X - gibbsDF$mu[i])^2)
  gibbsDF$tau[i] <- rgamma(1,shape = ag, rate = bg)
}

## @knitr resultsGibbs1

burnIn <- 1000

gibbsDFtidy <- gibbsDF %>%
  rowid_to_column(var = 'draw') %>%
  pivot_longer(-draw,values_to = 'value', names_to = 'parameter')


gibbsStat <- gibbsDFtidy %>%
  filter(draw >= burnIn) %>%
  group_by(parameter) %>%
  summarize(mean = mean(value),
            median = median(value),
            q2.5 = quantile(value,0.025),
            q97.5 = quantile(value,0.975))

mcmcPlot <- gibbsDFtidy %>%
  filter(draw >= burnIn) %>%
  ggplot(aes(x = draw, y = value, group = parameter, color = parameter)) + geom_line() + 
  facet_wrap(~parameter) + 
  labs(title = 'Trace Plot of Gibbs Sampler', x = 'Iteration') +
  theme(legend.position = 'none')
  
postPlot <- gibbsDFtidy %>%
  filter(draw >= burnIn) %>%
  ggplot(aes(x = value, fill = parameter)) + geom_density() + facet_wrap(~parameter) + 
  labs(title = 'Simulated Posterior') + 
  theme(legend.position = 'none')


## @knitr uberGibbs

# read data 
df <- read_rds('data/uberRides.rds') %>%
  mutate(amount = abs(amount),
         logAmount = log(amount))

minTrips <- 3

dfStats <- df %>%
  group_by(userId) %>%
  summarize(n_i = n(),
            sumY = sum(logAmount)) %>%
  filter(n_i >= minTrips) %>%
  mutate(userIdF = factor(userId),
         userIndex = as.integer(userIdF)) %>%
  ungroup() %>%
  arrange(userIndex)

dfIncl <- df %>%
  filter(userId %in% dfStats$userId)

ssqAll <- sum(dfIncl$logAmount^2)
nRiders <- nrow(dfStats)
nTotal <- nrow(dfIncl)

# set priors 
mu0 <- 0.0
tau0 <- 0.2
a0 <- 1.5
b0 <- 0.5
c0 <- 1.5
d0 <- 0.5

# start sampler 
nIter <- 3000
gibbs1 <- data.frame(mu_alpha = array(0,nIter), 
                      tau_alpha = array(0,nIter),
                      tau = array(0,nIter))
gibbs2 <- array(0,c(nIter,nRiders))

gibbs1$tau[1] <- 1.0        # initialization 
gibbs1$tau_alpha[1] <- 1.0  # initialization 
gibbs1$mu_alpha[1] <- 2.3   # initialization 

for (i in 2:nIter){
  ## alpha
  mu_var <- 1.0/(gibbs1$tau[i-1]*dfStats$n_i + gibbs1$tau_alpha[i-1])
  mu_mean <- mu_var*(gibbs1$tau[i-1]*dfStats$sumY + 
                      gibbs1$tau_alpha[i-1]*gibbs1$mu_alpha[i-1])
  gibbs2[i,] <- mu_mean + sqrt(mu_var)*rnorm(nRiders)
  
  ## mu_alpha
  mu_var <- 1.0/(gibbs1$tau_alpha[i-1]*nRiders + tau0)
  mu_mean <- mu_var*(gibbs1$tau_alpha[i-1]*sum(gibbs2[i,]) + tau0*mu0)
  gibbs1$mu_alpha[i] <- rnorm(1,mean = mu_mean, sd = sqrt(mu_var))
  
  ## tau_alpha
  ag <- c0 + 0.5*nRiders
  bg <- d0 + 0.5*sum((gibbs2[i,] - gibbs1$mu_alpha[i])^2)
  gibbs1$tau_alpha[i] <- rgamma(1,shape = ag, rate = bg)
  
  ## tau
  ag <- c0 + 0.5*nTotal
  bg <- d0 + 0.5*( ssqAll + sum(dfStats$n_i*(gibbs2[i,]^2)) - 
                     2.0*sum(gibbs2[i,]*dfStats$sumY))
  gibbs1$tau[i] <- rgamma(1,shape = ag, rate = bg)
  
}


burnIn <- 1000

gibbsDFtidy <- gibbs1 %>%
  rowid_to_column(var = 'draw') %>%
  pivot_longer(-draw,values_to = 'value', names_to = 'parameter')


gibbsStat <- gibbsDFtidy %>%
  filter(draw >= burnIn) %>%
  group_by(parameter) %>%
  summarize(mean = mean(value),
            median = median(value),
            q2.5 = quantile(value,0.025),
            q97.5 = quantile(value,0.975))

## @knitr uberGibbsPlots

tracePlotUber <- gibbsDFtidy %>%
  filter(draw >= burnIn) %>%
  ggplot(aes(x = draw, y = value, group = parameter, color = parameter)) + geom_line() + 
  facet_wrap(~parameter, scales = 'free') + 
  labs(title = 'Trace Plot of Gibbs Sampler', x = 'Iteration') +
  theme(legend.position = 'none')

# compare to Stan from week 2
thePars <- read_rds('../week 2/data/normalModelDraws.rds')

stanMCMCStat <- data.frame(mu_alpha = thePars$mu,
                       tau_alpha = 1.0/(thePars$sigma^2),
                       tau = 1.0/(thePars$sigma_y^2)) %>%
  rowid_to_column(var = 'draw') %>%
  pivot_longer(-draw,values_to = 'value', names_to = 'parameter') %>%
  group_by(parameter) %>%
  summarize(mean = mean(value),
            median = median(value),
            q2.5 = quantile(value,0.025),
            q97.5 = quantile(value,0.975))


alphaCompPlot <- data.frame(Stan = colMeans(thePars$alpha),
           Gibbs = colMeans(gibbs2[burnIn:nIter,])) %>%
  ggplot(aes(x = Stan, y = Gibbs)) + geom_point(alpha=0.5) + 
  labs(title = 'Posterior Mean of Alpha',
       subtitle = 'Comparison of Gibbs and Stan')








