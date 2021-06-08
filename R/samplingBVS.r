#' samples from the posterior distribution of a Bayesian Variable Selection model using weighted Tempered Gibbs Sampling
#' 
#' Perform Bayesian variable selection in linear regression contexts 
#' using discrete spike and slab priors. Posterior sampling and calculation 
#' of marginal Posterior Inclusion Probabilities (PIPs) for explanatory 
#' variables is done using the weighted Tempered Gibbs Sampling algorithm 
#' of Zanella and Roberts (2019).
#' 
#' The evaluated linear regression model can be written as
#' \deqn{Y|\beta_\gamma, \gamma, \sigma^2 ~ N(X_\gamma\beta_\gamma,\sigma^2(I_n))}
#' \deqn{\beta_\gamma|\gamma, \sigma^2 ~ N(0,\sigma^2\Sigma_\gamma)}
#' \deqn{p(\sigma^2) \propto 1/\sigma^2}
#' \deqn{\gamma_i|h iid~ Bern(h)       i = 1,...,p}
#' where the posterior probability of interest is \eqn{p(\gamma|Y)}. 
#' 
#' The prior covariance matrix of the coefficients 
#' of the selected regressors is \eqn{\Sigma_\gamma = c(X_\gamma^TX_\gamma)}, i.e. the g-prior recommended by Zellner 
#' (1986).
#' 
#' The Rao-Blackwellized estimators provide a vector with inclusion probabilities for 
#' each of the regressors, \eqn{{p(\gamma_i=1|Y)}_{i=1}^{p}}.
#' 
#' \eqn{h} can be a fixed value or \eqn{h ~ Beta(a,b)}.
#' 
#' The sampling algorithm flips one of p binary values of \eqn{\gamma} by sampling \eqn{i} 
#' from \eqn{1,...,p} proportionally to \eqn{p_i(\gamma)=p(\gamma_i|\gamma_{-i},Y)^{-1}}
#' in the case of Tempered Gibbs Sampling and proportionally to 
#' \eqn{p_i(\gamma)=(p(\gamma_i=1|\gamma_{-i},Y)+k/p)/p(\gamma_i|\gamma_{-i},Y)}
#' in the case of weighted Tempered Gibbs Sampling. Also, the weight of the new state of
#' the Markov Chain is proportional to \eqn{(\sum_{i=1}^p p_i(\gamma))^{-1}}.
#' 
#' For more information on weighted Tempered Gibbs Sampling, please refer to Zanella and
#' Roberts (2019).
#' 
#' 
#' @param y a vector of n observations (dependent variable) with dimensions (nx1).
#' @param X a matrix of p regressors (independent variables) with dimension (nxp).
#' @param c a real number greater than 0 which serves as a constant of proportionality
#' to the specification of the prior covariance matrix of the coefficients of the
#' selected regressors in the linear regression. The default is \code{NULL} which yields the recommended constant of proportionality for Zellner's g-prior 
#' , i.e. c = n.
#' @param h either a real number greater than 0 and smaller than 1 or a vector of real
#' values, both greater than 0. This parameter specifies the prior information of the
#' inclusion probability of the regressors which is identical for all regressors. 
#' In the former case, the prior probability is set to a fixed value. In the latter 
#' case, the prior probability is a Beta distribution with the specified parameters.
#' The default is the uniform distribution in terms of a Beta prior \code{c(1,1)}. 
#' @param n_iter a positive integer specifying the number of iterations for the Markov
#' Chain. The default is 2000.
#' @param burn_in either an integer greater than 1 or a real number greater than 0 and
#' smaller than 1. Specifies the number of burn in iterations for the Markov Chain. In
#' the former case the burn in iterations are set the fixed integer. In the latter case
#' the number of iterations are the specified percentage of the number of iterations.
#' The default is 0.2.
#' @param k_weight a real number greater than 0 which, in the case of \code{weighted = TRUE},
#' controls the tradeoff between exploration and exploitation in the choice of the variable
#' to be flipped at each iteration. A larger \code{k_weight} favours exploration. The default is 0.
#' @param weighted logical, with default \code{TRUE}, indicating whether to perform
#' weighted Tempered Gibbs Sampling if \code{TRUE} or Tempered Gibbs Sampling if 
#' \code{FALSE}.
#' 
#' @return A list with named objects:
#' \item{PIP }{a vector (px1) containing Rao-Blackwellised estimators of 
#' the marginal PIPs for each of the p regressors in \code{X}.}
#' \item{states }{a list containing the elements necessary to reproduce the samples of
#' the Markov Chain. These elements are:\cr
#' "start" - starting value for \eqn{\gamma} after the burnin period.\cr
#' "sample_weights" - a vector (n_iterx1) of weights for \eqn{\gamma} at each step of 
#' the Markov Chain.\cr
#' "indices_sequence" - a vector (n_iterx1) of indices ranging \eqn{1,...,p} indicating the
#' element of \eqn{\gamma} flipped at each step of the Markov Chain.}\cr
#' "y" - input independent variables.\cr
#' "X" - input dependent variables.\cr
#' "c" - constant of proportionality of the prior covariance matrix.
#'  
#' @export
#' 
#' @references
#' Zanella, G. and Roberts, G. (2019). Scalable importance tempering and Bayesian variable selection. Journal of the Royal Statistical Society: Series B (Statistical Methodology): 489–517. Crossref. Web.
#' 
#' Zellner, A. (1986). On Assessing Prior Distributions and Bayesian Regression Analysis with g-Prior Distributions. In: Goel, P. and Zellner, A., Eds., Bayesian Inference and Decision Techniques: Essays in Honor of Bruno de Finetti, Elsevier Science Publishers, Inc., New York, 233-243.
#' 
#' @seealso \code{\link{createSamples}} for creating the samples of the Markov Chain and their weights used to calculate the PIPs.
#' 
#' @examples
#' #Posterior inclusion probabilities of characteristics of cars on describing mileage
#' 
#' #load data
#' data(mtcars)
#' 
#' #create X matrix and y vector with zero mean for all regressors
#' X <- t(t(mtcars[,-1]) - colMeans(mtcars[,-1]))
#' y <- mtcars$mpg - mean(mtcars$mpg)
#' 
#' mtcars.output <- samplingBVS(y, X)
#' 
#' names(mtcars.output$PIP) <- names(mtcars[,-1])
#' print(mtcars.output$PIP)
samplingBVS <- function(y, #vector of observations
                        X, #matrix of regressors
                        c = NULL, #covariance matric and constant
                        h = c(1,1), #if vector 2x1 then parameters of Beta
                        n_iter = 2000, #number of effective iterations
                        burn_in = 0.2, #percentage(>0,<1)/number(>1) of burnin iterations
                        k_weight = 0, weighted = TRUE) { #weightedTGS and parameter
  
  ### wTGS algorithm for Bayesian variable selection problems
  
  ## Set options
  if (burn_in < 1) burn_in <- n_iter*burn_in
  if (any(h < 0)) stop("h must be a vector of two positive parameters of a Beta distribution or a real number between 0 and 1")
  
  ## throw errors for parameters
  if (burn_in < 0) stop("Burn in must be a real number between 0 and 1 or an integer larger than 1")
  if (!is.null(c)) 
    if (c <= 0) stop("c must be larger than 0")
  if (length(h) > 2 | length(h) == 0)
    stop("h must be a vector of two positive parameters of a Beta distribution or a real number between 0 and 1")
  else if (length(h) == 1)
    if (h > 1)
      stop("h must be a vector of two positive parameters of a Beta distribution or a real number between 0 and 1")

  if (n_iter < 0) stop("n_iter must be an positive integer")
  if (k_weight < 0) stop("k_weight must be larger than 0")
  
  ## check dimensions of y, X and the initial gamma
  n <- length(y)
  if (nrow(X) != n) stop("y and X should have the same number of observations")
  if (is.null(c)) c <- n 
  p <- ncol(X)
  
  if (length(h) == 2) {
	h1 <- h[1]
	h2 <- h[2]
  } else {
	h1 <- h
	h2 <- 0
  }
  
  output <- wTGS(as.matrix(X), as.vector(y), n, p, n_iter, burn_in, h1, h2, c, k_weight, weighted)
  
  return(list(PIP = output[[1]],
              states = list(start = output[[2]],
                            sample_weights = output[[3]],
                            indices_sequence = output[[4]]),
              y = as.vector(y),
              X = as.matrix(X),
              c = c))
}
