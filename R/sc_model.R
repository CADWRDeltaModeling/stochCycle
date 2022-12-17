

#' Build the stochastic cycle model plus trend matrices
#'
#' Constructs the matrices \code{FF} and \code{GG} used by DLM to represent
#' the stochastic cycle plus trend model
#' @param order_trend Order of the trend (low pass) component
#' @param order_cycle Order of the cyclical component, same for all frequencies
#' @param freqs frequencies deployed
#' @param rho dampening factor
#' @param trend_scale scaling factor trend to observation
#' @return list with attributes FF and GG
#' @export
sc_model_mat <- function(order_trend,order_cycle,freqs,rho,trend_scale=1.)
{
  # this function used to build the model
  freq_used <- freqs
  nfreq <- length(freq_used)
  nstate <- order_trend + 2*order_cycle*nfreq

  print(paste("number of frequencies=",nfreq))
  print(paste("order of trend=",order_trend))
  print(paste("order of cycle=",order_cycle))
  print(paste('number of state=',nstate))

  # Construct observation matrix FF

  obs_trend <- c(trend_scale,rep(0,order_trend-1))
  obs_per_freq <- c(1,0,rep(0,2*(order_cycle-1)))
  FF <- t(as.matrix(c(obs_trend, rep(obs_per_freq,nfreq))))


  # Construct state transition matrix GG
  # trend
  state_trend <- diag(order_trend)
  if (order_trend>1) {
    for (itrend in 1:(order_trend-1)){state_trend[itrend,itrend+1] = 1}
  }

  # cycles
  print("Stochastic cycle frequencies used")
  print(freq_used)
  print("Initial dampening (rho):")
  print(rho)

  single_freq <- function(f,rho){
    cycle_trig_part <- rep(list(rho*matrix(c(cos(f),-sin(f),sin(f),cos(f)),2)),order_cycle)
    single <- dlm::bdiag(cycle_trig_part)
    if (order_cycle > 1) {
      for(i in 1:(2*(order_cycle-1))){
        single[i,i+2] = 1
      }
    }
    single
  }

  cycle_matrices <-mapply(single_freq,freq_used,rho,SIMPLIFY=FALSE)

  state_cycle <- dlm::bdiag(cycle_matrices)

  # combine trend and cycle components
  GG <- dlm::bdiag(state_trend,state_cycle)
  #print("GG")
  #print (GG)
  list(FF=FF,GG=GG)
}


#' Return the last index of a stochastic cycle
#' @param order_trend order of the trend part of the model
#' @param order_cycle order of cycle part of model, assumed same for all freqs
#' @param nth_freq index of frequency being queried
find_cycle_index <- function(order_trend, order_cycle, nth_freq){
  phi1b_idx <- order_trend + nth_freq*order_cycle*2
  return (phi1b_idx)
}

#' Calculate sum of squared residuals from the state equations
#' corresponding to a single cycle according to equation (21) in
#' Harvey's paper (2005)
#' @param order_trend order of the trend part of the model
#' @param order_cycle order of cycle part of model, assumed same for all freqs
#' @param nth_freq index of frequency being queried
#' @param lambda frequency
#' @param rho dampening factor
#' @param alpha_draw state vector, for instance from posterior draw in Gibbs
#' @return ct the calculated sum squared residuals
#' @export
calc_ct_sum <- function(order_trend, order_cycle, nth_freq, lambda, rho, n, alpha_draw){
  phi1b_idx <- find_cycle_index(order_trend,order_cycle,nth_freq)
  phi1a_idx <- phi1b_idx - 1
  cos_lambda <- cos(lambda)
  sin_lambda <- sin(lambda)
  ct <- 0.0
  for (ii in (2:n)){
    ct <- ct + (alpha_draw[ii,phi1a_idx] - rho*cos_lambda*alpha_draw[(ii-1),phi1a_idx]
                - rho*sin_lambda*alpha_draw[(ii-1),phi1b_idx] )**2
    + (alpha_draw[ii,phi1b_idx] + rho*sin_lambda*alpha_draw[(ii-1),phi1a_idx]
       - rho*cos_lambda*alpha_draw[(ii-1),phi1b_idx] )**2
  }
  return(ct)
}

#' Shape parameter for the posterior distribution of kappa
#' according to equation (23) in Harvey's paper (2005). This is the sums of two
#' components, one from the unconditional covariance and the other from the sum of
#' squared residuals from the state equation. Beware this may be not exactly what this
#' function does
#' @param order_trend order of the trend part of the model
#' @param order_cycle order of cycle part of model, assumed same for all freqs
#' @param nth_freq index of frequency being queried
#' @param lambda frequency
#' @param rho dampening factor
#' @param covar_matrix_inv inverse of the covariance matrix (exactly?)
#' @return shape_kappa the shape function
#' @export
cal_shape_kappa <- function(order_trend, order_cycle, nth_freq, lambda, rho,
                            covar_matrix_inv, n, alpha_draw) {
  ct_sum <- calc_ct_sum(order_trend, order_cycle, nth_freq, lambda, rho, n, alpha_draw)
  cycle_idx_end <- find_cycle_index(order_trend, order_cycle, nth_freq)
  cycle_idx_start <- cycle_idx_end - 2*order_cycle + 1
  shape_kappa <- t(matrix(alpha_draw[1, (cycle_idx_start:cycle_idx_end)], ncol = 1)) %*%
    covar_matrix_inv %*%
    matrix(alpha_draw[1, (cycle_idx_start:cycle_idx_end)], ncol = 1)
  #todo: verbose
  #print(paste("shape kappa=",shape_kappa[1,1]," ct_sum=", ct_sum))
  shape_kappa=shape_kappa+ct_sum
  return (shape_kappa[1,1])   # the size is already 1x1, this just converts to scalar
}

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
sc_model_build <- function(order_trend,order_cycle,freqs,rho,sigma2_eps,sigma2_zeta,sigma2_kappa)
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
  C0_cycle <-  dlm::bdiag(covar)
  C0_full <- dlm::bdiag(C0_trend,C0_cycle)


  model <- dlm::dlm(FF=scmat$FF,V=V,GG=scmat$GG,W=W,m0=m0_full,C0=C0_full)

  return(model)
}

