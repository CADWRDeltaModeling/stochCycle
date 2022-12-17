

#' Calculate the filter lengths empirically using a unit impulse
#' @param order_trend Order of the trend (low pass) component
#' @param order_cycle Order of the cyclical component, same for all frequencies
#' @param freqs frequencies deployed
#' @param rho dampening factor
#' @param sigma2_eps variance of observation
#' @param sigma2_zeta variance of trend
#' @param sigma2_kappa variance of cycle disturbances
#' @return list of approximate trend extraction filter lengths (2*sigma)
#' @export
sc_filter_len <- function(order_trend,order_cycle,freqs,rho,sigma2_eps,sigma2_zeta,sigma2_kappa)
{
  nstate <- 2*order_cycle*length(freqs) + order_trend
  sigma2_diffuse <- 1.e7

  print(freqs)
  # prepare data for dlm with trend and cycle (Harvey 2005 paper)
  n_freq <- length(freqs)
  scmat <- sc_model_mat(order_trend, order_cycle, freqs, rho)

  # V is observation error, 1x1
  V <- matrix(sigma2_eps)

  # W is state innovation variance,

  Wtrend <- c(rep(0.,order_trend-1),1.)*sigma2_zeta
  Wcycle = c()
  for (jj in (1:n_freq)) {
    Wcycle <- c(Wcycle,rep(0,2*(order_cycle-1)),sigma2_kappa[jj],sigma2_kappa[jj])
  }
  W <- diag(c(Wtrend,Wcycle))
  m0_full <- rep(0, nstate)


  C0_trend <- diag(1,order_trend)*sigma2_diffuse
  covar <- list()
  for (jj in (1:n_freq)) {
    covar[[jj]] <- calc_cycle_covariance_matrix(order_cycle,rho,freqs[jj])*sigma2_kappa[jj]
  }
  C0_cycle <-  bdiag(covar)
  C0_full <- bdiag(C0_trend,C0_cycle)


  model <- dlm(FF=scmat$FF,V=V,GG=scmat$GG,W=W,m0=m0_full,C0=C0_full)

  pseudo <- c(rep(0,2000),c(1),rep(0,2000))
  modFilt <- dlmFilter(pseudo, model)
  modSmooth <- dlmSmooth(modFilt)
  return(modSmooth$s)
}
