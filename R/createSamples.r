#' creates the sample matrix from the posterior distribution of a Bayesian Variable Selection model
#' 
#' 
#' @param samples a list of states from the output of \code{samplingBVS} containing the elements necessary to reproduce the samples of
#' the Markov Chain. These elements are: \cr
#' "start" - starting value for \eqn{\gamma} after the burnin period. \cr
#' "sample_weights" - a vector (n_iterx1) of weights for \eqn{\gamma} at each step of 
#' the Markov Chain. \cr
#' "indices_sequence" - a vector (n_iterx1) of indices ranging \eqn{1,...,p} indicating the
#' element of \eqn{\gamma} flipped at each step of the Markov Chain.\cr
#' "y" - input independent variables.\cr
#' "X" - input dependent variables.\cr
#' "c" - constant of proportionality of the prior covariance matrix.
#' @param thin an integer greater than or equal to 1 indicating the period when to save
#' samples and sample weights from the Markov Chain, i.e. every how many steps of the 
#' Markov Chain the samples and sample weights should be recorded. The default is 1.
#'
#' @return A list with named objects:
#' \item{samples }{a matrix (n_iterxp) of \eqn{\gamma} at each step of the Markov Chain.}
#' \item{weights }{a vector (n_iterx1) of weights for \eqn{\gamma} at each step of 
#' the Markov Chain.}
#' \item{betas }{a matrix (n_iterxp) of \eqn{\beta} at each step of the Markov Chain.}
#' 
#' @export
#'
#' @seealso \code{\link{samplingBVS}} for running weighted Tempered Gibbs Sampling and caluculating Posterior Inclusion Probabiltiies.
#'
#' @examples
#' #Samples of inclusion of characteristics of cars on describing mileage
#' #load data
#' data(mtcars)
#' 
#' #create X matrix and y vector with zero mean for all regressors
#' X <- t(t(mtcars[,-1]) - colMeans(mtcars[,-1]))
#' y <- mtcars$mpg - mean(mtcars$mpg)
#' 
#' mtcars.output <- samplingBVS(y, X)
#' mtcars.samples <- createSamples(mtcars.output$states)
#' 
#' #Samples
#' head(mtcars.samples$samples)
createSamples <- function(samples, #list outputted from the main function
						              thin = 1) { #Thinning of the samples
  
  if (thin < 1) stop("thin must be an integer greater than or equal to 1")
  
  check <- names(samples)
  if (!all(check %in% c("start", "sample_weights", "indices_sequence", "y", "X", "c"))) {
	  if (!all(check %in% c("PIP", "states")))
	    stop("samples must be from the output of main sampling function (samplingBVS).")
	  else
	    samples <- samples$states
  }
  
  n_iter <- length(samples$indices_sequence)
  p <- length(samples$start)
  if (n_iter%/%thin*p > 10e8) warning(paste("Samples will be a", paste(n_iter%/%thin, p, sep = "x"), "matrix, consider thinning."))
  
  states <- matrix(NA, ncol = p, nrow = n_iter%/%thin)
  weights <- rep(NA, n_iter%/%thin)
  betas <- matrix(NA, ncol = p, nrow = n_iter%/%thin)
  
  gamma <- samples$start
  X <- samples$X
  XtX <- t(X) %*% X
  y <- samples$y
  n <- length(y)
  Xty <- c( t(X) %*% y )
  c <- samples$c
  for (t in 1:n_iter) {
    
    gamma[samples$indices_sequence[t]] <- 1 - gamma[samples$indices_sequence[t]]
    inv <- solve(XtX[gamma, gamma])
    b_hat <- c( inv %*% Xty[gamma] )
    s2 <- c( t(y - X %*% b_hat) %*% (y - X %*% b_hat) )
    sigma2 <- 1. / rgamma(1, n/2., s2/2. + 1./(2*(c+1)) * c( b_hat %*% Xty[gamma] ))
    beta <- MASS::mvrnorm(1, c/(c+1) * b_hat, sigma2 * c/(c+1) * inv)
    
    if (t %% thin == 0) {
      states[(t/thin),] <- gamma
      weights[t/thin] <- samples$sample_weights[t]
      beta_ <- rep(0, p)
      beta_[gamma] <- beta
      betas[(t/thin),] <- beta_
    }
  }
  
  return(list(samples = states,
              weights = weights,
              betas = betas))
}