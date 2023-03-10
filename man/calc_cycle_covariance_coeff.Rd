% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sc_uncond_cov.R
\name{calc_cycle_covariance_coeff}
\alias{calc_cycle_covariance_coeff}
\title{Calculate ij components for the unconditional covariance matrix
The coefficients are produced according to equation (8) in Trimbur's paper (2003) on higher orders stochastic cycles}
\usage{
calc_cycle_covariance_coeff(i, j, rho)
}
\arguments{
\item{i}{row of matrix}

\item{j}{column of matrix}

\item{rho}{dampening factor}
}
\value{
coefficient
}
\description{
Calculate ij components for the unconditional covariance matrix
The coefficients are produced according to equation (8) in Trimbur's paper (2003) on higher orders stochastic cycles
}
