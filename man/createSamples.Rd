% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createSamples.r
\name{createSamples}
\alias{createSamples}
\title{creates the sample matrix from the posterior distribution of a Bayesian Variable Selection model}
\usage{
createSamples(samples, thin = 1)
}
\arguments{
\item{samples}{a list of states from the output of \code{samplingBVS} containing the elements necessary to reproduce the samples of
the Markov Chain. These elements are: \cr
"start" - starting value for \eqn{\gamma} after the burnin period. \cr
"sample_weights" - a vector (n_iterx1) of weights for \eqn{\gamma} at each step of 
the Markov Chain. \cr
"indices_sequence" - a vector (n_iterx1) of indices ranging \eqn{1,...,p} indicating the
element of \eqn{\gamma} flipped at each step of the Markov Chain.}

\item{thin}{an integer greater than or equal to 1 indicating the period when to save
samples and sample weights from the Markov Chain, i.e. every how many steps of the 
Markov Chain the samples and sample weights should be recorded. The default is 1.}
}
\value{
A list with named objects:
\item{samples }{a matrix (n_iterxp) of \eqn{\gamma} at each step of the Markov Chain.}
\item{weights }{a vector (n_iterx1) of weights for \eqn{\gamma} at each step of 
the Markov Chain.}
}
\description{
creates the sample matrix from the posterior distribution of a Bayesian Variable Selection model
}
\examples{
#Samples of inclusion of characteristics of cars on describing mileage
#load data
data(mtcars)

#create X matrix and y vector with zero mean for all regressors
X <- t(t(mtcars[,-1]) - colMeans(mtcars[,-1]))
y <- mtcars$mpg - mean(mtcars$mpg)

mtcars.output <- samplingBVS(y, X)
mtcars.samples <- createSamples(mtcars.output$states)

#Samples
head(mtcars.samples$samples)
}
\seealso{
\code{\link{samplingBVS}} for running weighted Tempered Gibbs Sampling and caluculating Posterior Inclusion Probabiltiies.
}
