## @knitr loadLibraries
library(tidyverse)
library(tidybayes)
library(rstan)
library(ggridges)


## @knitr loadUberRides

df <- read_csv('data/uberRides.csv') %>%
  mutate(amount = abs(amount),
         logAmount = log(amount))

minTrips <- 3

dfStats <- df %>%
  group_by(userId) %>%
  summarize(nTrips = n(),
            avgLogPrice = mean(logAmount)) %>%
  filter(nTrips >= minTrips) %>%
  mutate(userIdF = factor(userId),
         userIndex = as.integer(userIdF)) %>%
  ungroup() %>%
  mutate(empiricalRank = n() + 1 - min_rank(avgLogPrice))

dfIncl <- df %>%
  filter(userId %in% dfStats$userId) %>%
  left_join(select(dfStats,userId,userIndex), by = 'userId')


## @knitr findTopRiders

extremes <- bind_rows(
  dfStats %>%
  slice_max(avgLogPrice,n=20) %>%
  mutate(Position = 'Top 20'),
  dfStats %>%
  slice_min(avgLogPrice,n=20) %>%
  mutate(Position = 'Bottom 20')
)


## @knitr plotUberStats

overallMean <- mean(dfStats$avgLogPrice)

dfStats %>%
  left_join(select(extremes,userId,Position), by = 'userId') %>%
  mutate(Position = replace_na(Position, 'Middle')) %>%
  ggplot(aes(x=nTrips,
             y=avgLogPrice,
             color = Position)) + geom_point(alpha=0.3) + geom_hline(aes(yintercept = overallMean)) +
  labs(x = 'Number of Trips',
       y = 'Average log(Trip Price)',
       title = 'Average Trip Price and Number of Trips') 


## @knitr uberModel

##
## Note: This chunk compiles models written in the probabilistic language Stan.
##       We will cover the details in week 4. For now we will just use the output  
##

options(mc.cores = parallel::detectCores())
normalModel <- stan_model(file='code/MLM_simple.stan')
saveRDS(normalModel, file = 'data/normalModel.rds')

normalModel <- read_rds('data/normalModel.rds')

Niter <-  1500
Nwarmup <- 1000
Nchains <- 4

fit <- sampling(normalModel,
                data=list(nObs = nrow(dfIncl),
                          nUsers = nrow(dfStats),
                          userID = dfIncl$userIndex,
                          y = dfIncl$logAmount),
                iter = Niter,
                warmup = Nwarmup,
                chains = Nchains,
                verbose = TRUE)
## save draws
thePars <- rstan::extract(fit,pars = c('alpha','mu','sigma','sigma_y'))
saveRDS(thePars, file = 'data/normalModelDraws.rds')


## @knitr getBayesEst1

thePars <- read_rds('data/normalModelDraws.rds')

## calculate posterior summary for each alpha and join in user info

# first get tidy version of all alpha draws 
alphaDraws <- as.data.frame(thePars$alpha) %>%
  rowid_to_column(var = 'draw') %>%
  pivot_longer(-draw,names_to = 'col',values_to = 'value') %>%
  separate(col, into = c('tmp','userIndex'),sep=1) %>%
  select(-tmp) %>%
  mutate(userIndex = as.integer(userIndex))


# calculate statistics and join user info
alphaStats <- alphaDraws %>%
  group_by(userIndex) %>%
  median_qi(value) %>%
  left_join(select(dfStats,userIndex,avgLogPrice,nTrips), by = 'userIndex')


## @knitr plotBayesEst1

alphaStats %>%
  select(value,avgLogPrice,nTrips) %>%
  rename(Bayes = value, Regular = avgLogPrice) %>%
  pivot_longer(names_to = 'Estimate',values_to = 'Value',-nTrips) %>%
  ggplot(aes(x=nTrips,
             y=Value,
             color = Estimate)
         ) + geom_point(alpha=0.3) + geom_hline(aes(yintercept = overallMean)) + facet_wrap(~Estimate) + 
  labs(x = 'Number of Trips',
       y = 'alpha',
       title = 'Bayes User Estimates and Number of Trips') 


## @knitr plotBayesEst2

alphaStats %>%
    select(value,avgLogPrice,nTrips) %>%
    ggplot(aes(x = avgLogPrice, y = value)) + geom_point(alpha=0.3, color ='red') + geom_abline(intercept = 0,slope = 1) + 
  labs(x = 'Regular Estimate',
       y = 'Bayes Estimate',
       title = 'Shrinkage Effect') 


## @knitr posteriorAlpha

# plot posterior of alpha for a few select users 

theUsers <- c(1143,1567,136,328)

alphaDraws %>%
  filter(userIndex %in% theUsers) %>%
  mutate(userIndex = factor(userIndex)) %>%
  ggplot(aes(x = value, fill = userIndex)) + geom_density(alpha=0.3) + 
  labs( x = 'alpha', title = 'Posterior for alpha for four users')


theStats <- alphaStats %>%
  filter(userIndex %in% theUsers) %>%
  select(-.point,-.interval,-.width)
  
  
## @knitr posteriorRanking

nUsers <- nrow(dfStats)

postRanks <- alphaDraws %>%
  ungroup() %>%
  group_by(draw) %>%
  mutate(rank = nUsers + 1 - min_rank(value))

postRanksStats <- postRanks %>%
  group_by(userIndex) %>%
  summarize(meanPostRank = mean(rank),
            minPostRank = min(rank),
            maxPostRank = max(rank)) %>%
  left_join(select(dfStats,userIndex,nTrips,avgLogPrice,empiricalRank), by = 'userIndex')


## @knitr retailData

retail <- read_csv('data/storeSales.csv') %>%
  mutate(storeIndex = as.integer(STORE_NUM))

storeID <- retail %>%
  select(STORE_NUM,storeIndex) %>%
  distinct(storeIndex,.keep_all = T)

overAllMeanLogPrice <- mean(retail$logFL)

retail$logFL_c <- retail$logFL - overAllMeanLogPrice

storeInfo <- read_csv('data/storeInfo.csv')

nObs <- nrow(retail)
nStores <- max(retail$storeIndex)

## @knitr StanModels

##
## Note: This chunk compiles models written in the probabilistic language Stan.
##       We will cover the details in week 4. For now we will just use the output  
##

## set options for processing markov chains in parallel
options(mc.cores = parallel::detectCores())

## compile and save Stan models
hiRegModel <- stan_model(file='code/MLM_regression.stan')
saveRDS(hiRegModel, file = 'data/hiRegModel.rds')

RegModelNoShrink <- stan_model(file='code/regression_no_shrink.stan')
saveRDS(RegModelNoShrink, file = 'data/RegModelNoShrink.rds')

hiRegModel2 <- stan_model(file='code/MLM_regression2.stan')
saveRDS(hiRegModel2, file = 'data/hiRegModel2.rds')



## @knitr RunSamplers 

## run Stan model with shrinkage 

hiRegModel <- read_rds('data/hiRegModel.rds')

Niter <-  1500
Nwarmup <- 1000
Nchains <- 4
 
fit <- sampling(hiRegModel,
               data=list(nObs = nObs,
                         nStores = nStores,
                         storeID = retail$storeIndex,
                         lp = retail$logFL_c,
                         ls = retail$logVol),
               iter = Niter,
               warmup = Nwarmup,
               chains = Nchains,
               verbose = TRUE)
## save draws
thePars <- rstan::extract(fit,pars = c('theta','mu','Omega','tau','sigma'))
saveRDS(thePars, file = 'data/hiRegModelDraws.rds')


## run Stan model with no shrinkage 
RegModelNoShrink <- read_rds('data/RegModelNoShrink.rds')

fitNoShrink <- sampling(RegModelNoShrink,
                        data=list(nObs = nObs,
                                  nStores = nStores,
                                  storeID = retail$storeIndex,
                                  lp = retail$logFL_c,
                                  ls = retail$logVol),
                        iter = Niter,
                        warmup = Nwarmup,
                        chains = Nchains,
                        verbose = TRUE)

## save draws
thePars <- rstan::extract(fitNoShrink,
                          pars = c('theta','sigma'))
saveRDS(thePars, file = 'data/RegModelNoShrinkDraws.rds')


## run Stan model with shrinkage and store level predictors 

zDF <- storeInfo %>%
  select(STORE_NUM,SALES_AREA_SIZE_NUM,SEG_VALUE_NAME) %>%
  mutate(intercept = 1) %>%
  mutate(STORE_NUM = factor(STORE_NUM),
         upscale = as.integer(SEG_VALUE_NAME == 'UPSCALE'),
         value = as.integer(SEG_VALUE_NAME == 'VALUE')) %>%
  left_join(storeID, by = 'STORE_NUM') %>%
  mutate(size = (SALES_AREA_SIZE_NUM - mean(SALES_AREA_SIZE_NUM))/1000) %>%
  arrange(storeIndex)

saveRDS(zDF, file = 'data/zDF.rds')

z <- zDF %>%
  select(intercept,size,upscale,value)

hiRegModel2 <- read_rds('data/hiRegModel2.rds')

Niter <-  1500
Nwarmup <- 1000
Nchains <- 4

fit <- sampling(hiRegModel2,
                data=list(nObs = nObs,
                          nStores = nStores,
                          dimZ = ncol(z),
                          storeID = retail$storeIndex,
                          lp = retail$logFL_c,
                          ls = retail$logVol,
                          Z = as.matrix(z)),
                iter = Niter,
                warmup = Nwarmup,
                chains = Nchains,
                verbose = TRUE)
## save draws
thePars <- rstan::extract(fit,pars = c('theta','gamma','Omega','tau','sigma'))

saveRDS(thePars, file = 'data/hiRegModel2Draws.rds')



## @knitr retailResultsShrinkage

## read draws 
thePars <- read_rds('data/hiRegModelDraws.rds')

## convert draws to tidy format 
alphaDraws <- as.data.frame(thePars$theta[,,1]) %>%
  rowid_to_column(var = 'draw') %>%
  pivot_longer(-draw,names_to = 'var',values_to = 'value') %>%
  separate(var,into = c('tmp1','storeIndex'),sep=1) %>%
  select(-tmp1)

betaDraws <- as.data.frame(thePars$theta[,,2]) %>%
  rowid_to_column(var = 'draw') %>%
  pivot_longer(-draw,names_to = 'var',values_to = 'value') %>%
  separate(var,into = c('tmp1','storeIndex'),sep=1) %>%
  select(-tmp1)

## calculate posterior statistics for slopes and intercepts 
betaStats <- betaDraws %>%
  group_by(storeIndex) %>%
  mean_qi(value)

alphaStats <- alphaDraws %>%
  group_by(storeIndex) %>%
  mean_qi(value)


## @knitr plotResultsShrinkage

alphaStatsPlot <- alphaStats %>%
  mutate(storeIndex = factor(storeIndex)) %>%
  ggplot(aes(x = fct_reorder(storeIndex,value), y = value, ymin = .lower, ymax = .upper)) + 
  geom_pointrange() + 
  labs(x = 'Store Index', 
       y = expression(alpha),
       subtitle = 'Posterior Mean and 95% inner quantile range',
       title = expression(paste('Posterior Summary for ', alpha))) + 
  theme(axis.text.x  = element_text(size=6))

betaStatsPlot <- betaStats %>%
  mutate(storeIndex = factor(storeIndex)) %>%
  ggplot(aes(x = fct_reorder(storeIndex,value), y = value, ymin = .lower, ymax = .upper)) + 
  geom_pointrange() + 
  labs(x = 'Store Index', 
       y = expression(beta),
       subtitle = 'Posterior Mean and 95% inner quantile range',
       title = expression(paste('Posterior Summary for ', beta))) +
  theme(axis.text.x  = element_text(size=6))



alphaStatsOrder <- alphaStats %>%
  arrange(value)

alphaFullPosterior <- alphaDraws %>%
  mutate(storeIndex = factor(storeIndex, levels = alphaStatsOrder$storeIndex)) %>%
  ggplot(aes(x = value, y = storeIndex)) + 
  geom_density_ridges(alpha=0.5, rel_min_height = 0.01) + 
  labs(y = 'Store Index', 
       x = expression(alpha),
       title = expression(paste('Posterior for ', alpha))) +
  theme_bw() + 
  theme(axis.text.y  = element_text(size=6))

betaStatsOrder <- betaStats %>%
  arrange(value)

betaFullPosterior <- betaDraws %>%
  mutate(storeIndex = factor(storeIndex, levels = betaStatsOrder$storeIndex)) %>%
  ggplot(aes(x = value, y = storeIndex)) + 
  geom_density_ridges(alpha=0.5, rel_min_height = 0.01) + 
  labs(y = 'Store Index', 
       x = expression(beta),
       title = expression(paste('Posterior for ', beta))) +
  theme_bw() + 
  theme(axis.text.y  = element_text(size=6))

  

jointStats <- data.frame(alpha = alphaStats$value,
                         beta = betaStats$value,
                         storeIndex = as.integer(betaStats$storeIndex)) %>%
  ggplot(aes(x = alpha, y = beta, label = storeIndex)) + geom_label() +
  geom_smooth(method='lm',se = FALSE) + 
  labs(title = expression(paste('Posterior Mean of ',alpha, ' and ', beta)),
       subtitle = 'Posterior Means for 76 Stores',
       x = expression(alpha),
       y = expression(beta))


theStores <- c("56","15","25","67","50")

theAlphas <- alphaDraws %>%
  filter(storeIndex %in% theStores) %>%
  rename(alpha = value)

theBetas <- betaDraws %>%
  filter(storeIndex %in% theStores) %>%
  rename(beta = value)

plotFull <- theAlphas %>%
  left_join(theBetas, by = c('draw','storeIndex')) %>%
  ggplot(aes(x = alpha,y = beta, color = storeIndex)) + 
  geom_point(alpha=0.2,size=0.3) +
  guides(colour = guide_legend(override.aes = list(size=2,alpha=1))) + 
  labs(title = expression(paste('Joint Posterior of ',alpha, ' and ', beta)),
       subtitle = '2,000 Posterior Draws for each of 5 Stores',
       x = expression(alpha),
       y = expression(beta))


## @knitr plotResultsNoShrinkage

thePars <- read_rds('data/RegModelNoShrinkDraws.rds')

## convert draws to tidy format 
alphaDrawsNoShrink <- as.data.frame(thePars$theta[,,1]) %>%
  rowid_to_column(var = 'draw') %>%
  pivot_longer(-draw,names_to = 'var',values_to = 'value') %>%
  separate(var,into = c('tmp1','storeIndex'),sep=1) %>%
  select(-tmp1)

betaDrawsNoShrink <- as.data.frame(thePars$theta[,,2]) %>%
  rowid_to_column(var = 'draw') %>%
  pivot_longer(-draw,names_to = 'var',values_to = 'value') %>%
  separate(var,into = c('tmp1','storeIndex'),sep=1) %>%
  select(-tmp1)

betaStatsNoShrink <- betaDrawsNoShrink %>%
  group_by(storeIndex) %>%
  mean_qi(value) %>%
  mutate(model = 'No Shrinkage')

alphaStatsNoShrink <- alphaDrawsNoShrink %>%
  group_by(storeIndex) %>%
  mean_qi(value) %>%
  mutate(model = 'No Shrinkage')

betaStatsBoth <- bind_rows(
  betaStats %>%
    mutate(model = 'Shrinkage'),
  betaStatsNoShrink
)


nStorePlot <- 30

storeOrder <- betaStats %>%
  sample_n(nStorePlot) %>%
  arrange(value) %>%
  select(storeIndex)

betaStatsBoth %>%
  filter(storeIndex %in% storeOrder$storeIndex) %>%
  mutate(storeIndex = factor(storeIndex, levels = storeOrder$storeIndex)) %>%
  ggplot(aes(x = model, y = value, ymin = .lower, ymax = .upper, color = model)) +
  geom_pointrange() + 
  facet_wrap(~storeIndex, nrow=1) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none") + 
  labs(title = 'Posterior Summaries with Shrinkage and No Shrinkage',
       subtitle = paste0(nStorePlot,' random stores'))







## @knitr plotExpectedDemand

theAlphas <- alphaDraws %>%
  rename(alpha = value)

theBetas <- betaDraws %>%
  rename(beta = value)

theSigmas <- data.frame(sigma = thePars$sigma) %>%
  rowid_to_column(var = 'draw')
  
allDraws <- theAlphas %>%
  left_join(theBetas, by = c('draw','storeIndex')) %>%
  left_join(theSigmas, by = c('draw'))

theStores <- c("56","15")

allDrawsTheStores <- allDraws %>%
  filter(storeIndex %in% theStores)


theBetaPoints <- betaStats %>% 
  filter(storeIndex %in% theStores)

theAlphaPoints <- alphaStats %>% 
  filter(storeIndex %in% theStores)

sigmaPoint <- mean(thePars$sigma)


priceSeq <- c(3,3.25,3.5,3.75)

demandDraws <- NULL
demandPoint <- NULL

for (price in priceSeq){

  theMeanDraws <- exp(allDrawsTheStores$alpha + allDrawsTheStores$beta*(log(price) - overAllMeanLogPrice) + 
                        0.5*(allDrawsTheStores$sigma^2))

  demandDraws <- bind_rows(demandDraws, data.frame(demand = theMeanDraws, 
                                                   price = price,
                                                   storeIndex = allDrawsTheStores$storeIndex))
  
  demandPoint <- bind_rows(demandPoint,
                           data.frame(
                             demandPoint = exp(theAlphaPoints$value + 
                                                 theBetaPoints$value*(log(price) - overAllMeanLogPrice) + 
                                                 0.5*(sigmaPoint^2)),
                             price = price,
                             storeIndex = theAlphaPoints$storeIndex)
                           )
  
}




tmp1 <- demandDraws %>%
  mutate(price = factor(price),
         storeLabel = paste0('Store Index ',storeIndex))

tmp2 <- demandPoint %>%
  mutate(price = factor(price),
         storeLabel = paste0('Store Index ',storeIndex))

ggplot(data = tmp1, aes(x = demand, fill = price)) + geom_density(alpha=0.4) + 
  geom_vline(data = tmp2, aes(xintercept = demandPoint, color = price),show.legend = F) + 
  facet_wrap(~storeLabel, nrow = 1) + 
  labs(title = 'Posterior of Expected Demand',
       subtitle = 'Two Stores',
       caption = 'Note: Vertical lines indicate expected demand at posterior mean values')
  


## @knitr posteriorPredictive 

theStore <- "15"
allDrawsTheStore <- allDrawsTheStores %>%
  filter(storeIndex == theStore)

demandDraws <- NULL

nDraws <- nrow(allDrawsTheStore)

for (price in priceSeq){
  
  theMeans <- allDrawsTheStore$alpha + allDrawsTheStore$beta*(log(price) - overAllMeanLogPrice)
  theSDs <- allDrawsTheStore$sigma
  theDraws <- rlnorm(nDraws, meanlog = theMeans, sdlog = theSDs)
  
  demandDraws <- bind_rows(demandDraws, data.frame(demand = theDraws, 
                                                   price = price)
  )

}

demandDraws %>%
  mutate(price = factor(price)) %>%
  ggplot(aes(x = demand, fill = price)) + geom_density(alpha=0.4) + 
  labs(title = 'Posterior Predictive Demand Distribution',
       subtitle = paste0('Store Index ',theStore))


## @knitr resultsModel2

thePars <- read_rds('data/hiRegModel2Draws.rds')
zDF <- read_csv('data/zDF.csv')




zNames <- c('Intercept','Store Size','Upscale?','Value?')
parsDF <- bind_rows(
  as.data.frame(thePars$gamma[,,1]) %>%
    mutate(parameter = 'alpha'),
  as.data.frame(thePars$gamma[,,2]) %>%
    mutate(parameter = 'beta')
)
  
names(parsDF)[1:4] <- zNames

parsTidy <- parsDF %>%
  pivot_longer(-parameter, names_to = 'Effect', values_to = 'Value')

alphaPlot <- parsTidy %>%
  filter(parameter == 'alpha') %>%
  ggplot(aes(x = Value, fill = Effect)) + geom_density(alpha=0.5) + 
  facet_wrap(~Effect, scales = 'free', nrow = 2) + 
  labs(title = 'Gamma Posterior for alpha equation')

betaPlot <- parsTidy %>%
  filter(parameter == 'beta') %>%
  ggplot(aes(x = Value, fill = Effect)) + geom_density(alpha=0.5) + 
  facet_wrap(~Effect, scales = 'free', nrow = 2) + 
  labs(title = 'Gamma Posterior for beta equation')


## scatter plot of effects with store information
alphaStats <- as.data.frame(thePars$theta[,,1]) %>%
  rowid_to_column(var = 'draw') %>%
  pivot_longer(-draw,names_to = 'var',values_to = 'value') %>%
  separate(var,into = c('tmp1','storeIndex'),sep=1) %>%
  select(-tmp1) %>%
  group_by(storeIndex) %>%
  mean_qi(value) %>%
  rename(alpha = value) %>%
  select(storeIndex,alpha)

allStats <- as.data.frame(thePars$theta[,,2]) %>%
  rowid_to_column(var = 'draw') %>%
  pivot_longer(-draw,names_to = 'var',values_to = 'value') %>%
  separate(var,into = c('tmp1','storeIndex'),sep=1) %>%
  select(-tmp1) %>%
  group_by(storeIndex) %>%
  mean_qi(value) %>%
  rename(beta = value) %>%
  select(storeIndex,beta) %>%
  left_join(alphaStats, by = 'storeIndex') %>%
  mutate(storeIndex = as.integer(storeIndex)) %>%
  left_join(zDF, by = 'storeIndex')


alphaBetaPlot <- allStats %>%
  ggplot(aes(x = alpha, y = beta, color = SEG_VALUE_NAME)) + geom_point() + 
  labs(title = 'Posterior Estimates and Store Segment')





