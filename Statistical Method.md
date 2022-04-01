# Statistical Method

Created: March 28, 2022 6:32 PM
Last Edited Time: March 28, 2022 6:57 PM
Status: In Progress

## Steps of Bayesian data analysis

1. Setting up a *full probability model*—a joint probability distribution for all observable and unobservable quantities in a problem. The model should be consistent with knowledge about the underlying scientific problem and the data collection process.
2. Conditioning on observed data: calculating and interpreting the appropriate posterior distribution—the conditional probability distribution of the unobserved quantities of ultimate interest, given the observed data.
3. Evaluating the fit of the model and the implications of the resulting posterior distribution: how well does the model fit the data, are the substantive conclusions reasonable, and how sensitive are the results to the modeling assumptions in step 1? In response, one can alter or expand the model and repeat the three steps.

## Parameters, notations

$\theta$  is unobservable vector quantities or population parameters of interest (such as the probabilities of survival under each treatment for randomly chosen members of the population in the example of the clinical trial)

$y$ denote the observed data (such as the numbers of survivors and deaths in each treatment group)

$\tilde{y}$ is unknown but potentially observable

$u^Tu$ is a scalar and $uu^T$ is an n×n matrix

## Bayes’ Rule

先验概率是不加信息判断一件事的概率，比如是否聪明的概率，是否患病的概率。后验概率是加入一定信息的概率，是条件概率。比如我知道你考上好大学，那么理应比平均智商高些。

$p(\theta, y)=p(\theta)p(y|\theta)$

$p(\theta|y) = \frac{p(\theta, y)}{p(y)} = \frac{p(\theta)p(y|\theta)}{p(y)}$, where $p(y) = \sum_{\theta}p(\theta)p(y|\theta)$


e.g. $Pr(\theta>2) = \int_{\theta>2}p(\theta)d\theta$

variance is defined as $sd(\theta)/E(\theta)$, geometric mean is $exp(E[log(\theta)])$

