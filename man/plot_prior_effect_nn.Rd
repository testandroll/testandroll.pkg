% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nn_functions.R
\name{plot_prior_effect_nn}
\alias{plot_prior_effect_nn}
\title{Plots prior densities against treatment effect (difference in profit per customer)}
\usage{
plot_prior_effect_nn(mu, sigma, abs = FALSE)
}
\arguments{
\item{mu}{means of the prior on the mean response}

\item{sigma}{standard deviations of the prior on the mean response}

\item{abs}{whether or not to take the absolute difference, 'TRUE' or 'FALSE'}
}
\value{
graph plotting prior density vs. difference in profit per customer
}
\description{
Plots prior densities against treatment effect (difference in profit per customer)
}
\examples{
plot_prior_effect_nn(mu=c(5,10), sigma=c(10,10), abs=FALSE)
}
