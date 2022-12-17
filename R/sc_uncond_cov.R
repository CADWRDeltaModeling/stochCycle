#' Calculate ij components for the unconditional covariance matrix
#' The coefficients are produced according to equation (8) in Trimbur's paper (2003) on higher orders stochastic cycles
#'
#' @param i row of matrix
#' @param j column of matrix
#' @param rho dampening factor
#' @return coefficient
calc_cycle_covariance_coeff <- function(i,j,rho) {
  if (i > j) {
    temp <- j
    j <- i
    i <- temp }
  denominator <- (1-rho**2)**(i+j-1)
  i_minus_1 <- i -1
  numerator <- 0
  for (r in (0:i_minus_1)) {
    numerator <- numerator + choose(i_minus_1,r)*choose((j-1),(r+j-i))*rho**(2*r)
  }
  return (numerator/denominator)
}

#' Calculate ij components for the unconditional covariance matrix
#' The coefficients are produced according to equation (8) in Trimbur (2003)
#' on higher orders stochastic cycles. What is difference with prior?
#'
#' @param i row of matrix
#' @param j column of matrix
#' @param order_cycle the order (n) of the stochastic cycle model
#' @param rho dampening factor
#' @param lambda frequency of the cycle
#' @return ij entry of matrix
calc_cycle_variance_matrix_ij <- function(i,j,order_cycle,rho,lambda){
  # this function calculates values of ij components for the covariance matrix
  # according to equation (8) in Trimbur's paper (2003) on higher orders stochastic cycles
  multiplier <- calc_cycle_covariance_coeff(i,j,rho)
  do_transpose <- FALSE
  temp1 <- rho*cos(lambda)
  temp2 <- rho*sin(lambda)
  T <- matrix(c(temp1,-1*temp2,temp2,temp1), nrow=2)
  if (i == j) {
    return ((multiplier*diag(1,2)))}
  else {
    if (i > j) {
      do_transpose <- TRUE
      temp <- j
      j <- i
      i <- temp }
    temp <- diag(1, 2)
    for (k in (1:(j - i))) {
      temp <- temp %*% T
    }
    if (do_transpose) {
      return(t(multiplier*temp))
    }
    else {
      return((multiplier*temp))
    }
  }
}

#' Build unconditional covariance matrix for cycle of one frequency
#' @param order_cycle order of cycle
#' @param rho dampening factor
#' @param lambda frequency of cycle
#' @return covariance matrix
#' @export
calc_cycle_covariance_matrix <- function(order_cycle,rho,lambda) {
  if (order_cycle == 1) {
    covariance_matrix <- calc_cycle_variance_matrix_ij(1,1,1,rho,lambda)
  }
  else{
    for (ii in (order_cycle:1)) {
      for (jj in (order_cycle:1)) {
        temp <- calc_cycle_variance_matrix_ij(ii, jj, order_cycle, rho, lambda)
        if (jj == order_cycle) {temp2 <- temp}
        else {temp2 <- cbind(temp2, temp)}
      }
      if (ii==order_cycle) {temp3 <- temp2 }
      else {temp3 <- rbind(temp3, temp2)}
    }
    covariance_matrix <- temp3
  }
  if (!lqmm::is.positive.definite(covariance_matrix,method='chol'))
  {print ("Warning: Covariance isn't positive definite!")}
  return(covariance_matrix)
}
