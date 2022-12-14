% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nn_functions.R
\name{test_size_nht}
\alias{test_size_nht}
\title{Computes the recommended sample sizes for a null hypothesis test}
\usage{
test_size_nht(s, d, conf = 0.95, power = 0.8, N = NULL)
}
\arguments{
\item{s}{vector of length 1 (symmetric) or 2 (asymmetric) indicating
response standard deviations}

\item{d}{minimum detectable difference between treatments}

\item{conf}{1 - type I error rate}

\item{power}{1 - type II error rate}

\item{N}{finite deployment population, if NULL no finite population correction is used}
}
\value{
recommended sample sizes as a vector
}
\description{
Computes the recommended sample sizes for a null hypothesis test
}
\examples{
test_size_nht(s=c(0.5,0.10), d=0.2, conf=0.95, power=0.8, N=NULL)
test_size_nht(s=0.5, d=0.2, conf=0.95, power=0.8, N=NULL)
}
