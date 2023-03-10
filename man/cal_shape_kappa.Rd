% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sc_model.R
\name{cal_shape_kappa}
\alias{cal_shape_kappa}
\title{Shape parameter for the posterior distribution of kappa
according to equation (23) in Harvey's paper (2005). This is the sums of two
components, one from the unconditional covariance and the other from the sum of
squared residuals from the state equation. Beware this may be not exactly what this
function does}
\usage{
cal_shape_kappa(
  order_trend,
  order_cycle,
  nth_freq,
  lambda,
  rho,
  covar_matrix_inv,
  n,
  alpha_draw
)
}
\arguments{
\item{order_trend}{order of the trend part of the model}

\item{order_cycle}{order of cycle part of model, assumed same for all freqs}

\item{nth_freq}{index of frequency being queried}

\item{lambda}{frequency}

\item{rho}{dampening factor}

\item{covar_matrix_inv}{inverse of the covariance matrix (exactly?)}
}
\value{
shape_kappa the shape function
}
\description{
Shape parameter for the posterior distribution of kappa
according to equation (23) in Harvey's paper (2005). This is the sums of two
components, one from the unconditional covariance and the other from the sum of
squared residuals from the state equation. Beware this may be not exactly what this
function does
}
