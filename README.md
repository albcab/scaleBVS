
<!-- README.md is generated from README.Rmd. Please edit that file -->

# [scaleBVS](https://github.com/albcab/scaleBVS)

Bayesian Variable Selection models provide a natural and coherent
framework to select a subset of explanatory variables in linear
regression contexts. The binary inclusion variables for each regressor
typically possess pairwise and/or negative dependence structures
conjectured to be conductive to successful application of weighted
Tempered Gibbs Sampling (Zanella and Roberts, 2019).

The use of weighted Tempered Gibbs Sampling overcomes the challenges of
high-dimensional Bayesian Variable selection models by an efficient
computation of the full conditional distribution of the binary inclusion
probabilities. These full conditionals allow for the calculation of
Rao-Blackwellised estimators of the marginal Posterior Inclusion
Probabilities for each regressor. These estimates quantify the
uncertainties of the true underlying linear model.

This package has been concieved as an implementation of the weighted
Tempered Gibbs Sampling algorithm of Zanella and Roberts (2019) to
Bayesian Variable Selection models in order to sample from the
distribution of its binary inclusion variables and provide a formal
Bayesian answer to variable selection problems.

## Installation

Install scaleBVS from [Github](https://github.com/albcab/scaleBVS) with:

``` r
install.packages("devtools")
devtools::install_github("albcab/scaleBVS")
```

You can also track or contribute to the development of scaleBVS at
<https://github.com/albcab/scaleBVS>.

## Features

scaleBVS samples from the posterior distribution of the vector of binary
inclusion variables
<img src="https://render.githubusercontent.com/render/math?math=\gamma = (\gamma_1,...,\gamma_p) \in \{1,0\}^p">
for the linear model with p regressors and n
observations:

<img src="https://render.githubusercontent.com/render/math?math=\mathbf{Y}|\beta_{\gamma},\gamma,\sigma^2 \sim N(X_{\gamma}\beta_{\gamma},\sigma I_n)">

<img src="https://render.githubusercontent.com/render/math?math=\beta_{\gamma}|\gamma,\sigma^2 \sim N(0,c(X_{\gamma}^TX_{\gamma})^{-1})">

<img src="https://render.githubusercontent.com/render/math?math=\gamma_i|h \overset{iid}{\sim} Bernoulli(h)">

where
<img src="https://render.githubusercontent.com/render/math?math=c > 0">
and
<img src="https://render.githubusercontent.com/render/math?math=\sigma^2">
is given a noninformative prior. The variable
<img src="https://render.githubusercontent.com/render/math?math=h"> can
be either set to a fixed value between 0 and 1 or given a Beta prior
distribution with two fixed hyperparameters, both greater than 0.

The algorithm has a uniform and a weighted version. The latter allows
for a tradeoff between exploration, forces sampler to explore new
regions, and exploitation, sampler focuses on most important regions, of
the state space.

Finally, in the case that
<img src="https://render.githubusercontent.com/render/math?math=X_{\gamma}^TX_{\gamma}">
is singular, the Moore-Penrose inverse (pseudo or generalized inverse)
of the matrix is used.

For more information on weighted Tempered Gibbs Sampling, please refer
to Zanella and Roberts (2019).

## Usage

scaleBVS contains one main function `samplingBVS()` used to calculate
Posterior Inclusion Probabilities for any set of p regressors. The
function `createSamples` computes the matrix of samples and weights
given the output of `samplingBVS`. Since the model does not include an
intercept on the regessors, it is recommended to standardize both the
dependent and independent variables to a zero mean.

``` r
#Standardizing
X <- t(t(mtcars[,-1]) - colMeans(mtcars[,-1]))
y <- mtcars$mpg - mean(mtcars$mpg)

#SAMPLING
mtcars.output <- scaleBVS::samplingBVS(y, X)

#POSTERIOR INCLUSION PROBABILITIES
names(mtcars.output$PIP) <- names(mtcars[,-1])
print(mtcars.output$PIP)
#>       cyl      disp        hp      drat        wt      qsec        vs 
#> 0.3723372 0.1544439 0.3671071 0.1446202 0.9224017 0.3535772 0.1336101 
#>        am      gear      carb 
#> 0.2595719 0.1402785 0.2077849

#SAMPLES
mtcars.samples <- scaleBVS::createSamples(mtcars.output$states)
head(mtcars.samples$samples, 3)
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> [1,]    0    0    0    1    1    0    0    0    0     0
#> [2,]    0    0    0    1    0    0    0    0    0     0
#> [3,]    0    0    0    1    1    0    0    0    0     0
```

## Benchmarks

scaleBVS is very fast and can easily handle high-dimensional problems.
We showcase these characteristics through an illustration using
simulated data. The simulated dataset has a very large amount of
uncorrelated regressors where only the first 10 are strongly correlated
with the response variable.

``` r
#Setting up the number of observations, paramters and samples
n <- 10000
p <- 10000
n_iter <- 1000
burnin <- 200

#Creating the data
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
beta0 <- matrix(c(2,-3,2,2,-3,3,-2,3,-2,3, rep(0, p-10)), ncol = 1)
beta <- 2*sqrt(log(p)/n)*beta0
y <- X %*% beta + rnorm(n, mean = 0, sd = 1)
X <- t(t(X) - colMeans(X))
y <- y - mean(y)

#TIMING
microbenchmark::microbenchmark(
  scaleBVS = scaleBVS::samplingBVS(y, X, n_iter = n_iter, burn_in = burnin),
  times = 10
)
#> Unit: seconds
#>      expr      min       lq     mean   median       uq      max neval
#>  scaleBVS 3.680549 4.404221 5.074466 4.569247 5.132157 9.638155    10

#RESULTS CHECK
head(scaleBVS::samplingBVS(y, X, n_iter = n_iter, burn_in = burnin)$PIP, 11)
#>  [1] 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00
#>  [6] 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00
#> [11] 2.254658e-05
```

The previous example is run on a Windows server 64-bit operating system
with 4GB of RAM and Intel i5 CPU. With more memory `samplingBVS()` can
easily handle even larger problems.

## References

  - Zanella, G. and Roberts, G. (2019). Scalable importance tempering
    and Bayesian variable selection. Journal of the Royal Statistical
    Society: Series B (Statistical Methodology): 489–517. Crossref. Web.
  - Zellner, A. (1986). On Assessing Prior Distributions and Bayesian
    Regression Analysis with g-Prior Distributions. In: Goel, P. and
    Zellner, A., Eds., Bayesian Inference and Decision Techniques:
    Essays in Honor of Bruno de Finetti, Elsevier Science Publishers,
    Inc., New York, 233-243.
