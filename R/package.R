#' Bayesian Variable Selection for the linear model via weighted Tempered Gibbs Sampling
#' 
#' Perform Bayesian variable selection in linear regression contexts using discrete 
#' spike and slab priors. Posterior sampling and calculation of marginal Posterior 
#' Inclusion Probabilities (PIPs) for explanatory variables is done using the 
#' weighted Tempered Gibbs Sampling algorithm of Zanella and Roberts (2019).
#' 
#' Bayesian Variable Selection models provide a natural and coherent framework
#' to select a subset of explanatory variables in linear regression contexts. 
#' The binary inclusion variables for each regressor typically possess pairwise
#' and/or negative dependence structures conjectured to be conductive to 
#' successful application of weighted Tempered Gibbs Sampling (Zanella and Roberts, 2019).
#' 
#' The use of weighted Tempered Gibbs Sampling overcomes the challenges of
#' high-dimensional Bayesian Variable selection models by an efficient computation 
#' of the full conditional distribution of the binary inclusion probabilities.
#' These full conditionals allow for the calculation of Rao-Blackwellised 
#' estimators of the marginal Posterior Inclusion Probabilities for each regressor.
#' These estimates quantify the uncertainties of the true underlying
#' linear model.
#' 
#' This package has been concieved as an implementation of the weighted
#' Tempered Gibbs Sampling algorithm to Bayesian Variable Selection models in order
#' to sample from the distribution of its binary inclusion variables and provide
#' a formal Bayesian answer to variable selection problems.
#' 
#' \tabular{ll}{ Package: \tab scaleBVS\cr Type: \tab Package\cr Version:
#' \tab 1.0.0\cr Date: \tab 2020-01-20\cr License: \tab GPL-2\cr }
#' 
#' @name scaleBVS-package
#' @aliases scaleBVS-package scaleBVS
#' @docType package
#' @author Giacomo Zanella and Alberto Cabezas Gonzalez
#' 
#' Maintainer: Alberto Cabezas Gonzalez \email{alb.cab94@gmail.com}
#' 
#' @seealso \code{\link{samplingBVS}}, \code{\link{createSamples}}
#' 
#' @references
#' Zanella, G. and Roberts, G. (2019). Scalable importance tempering and Bayesian variable selection. Journal of the Royal Statistical Society: Series B (Statistical Methodology): 489–517. Crossref. Web.
#' 
#' Zellner, A. (1986). On Assessing Prior Distributions and Bayesian Regression Analysis with g-Prior Distributions. In: Goel, P. and Zellner, A., Eds., Bayesian Inference and Decision Techniques: Essays in Honor of Bruno de Finetti, Elsevier Science Publishers, Inc., New York, 233-243.
#' 
#' @keywords package
NULL