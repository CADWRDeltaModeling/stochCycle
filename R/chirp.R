

Qr_orig <- function(x,c1=3.0,c2=5.5){
  3.5+5.5*AiryA(1.-.0002*(x-4.75*pi)^2.)
}


Qr_reflect <- function(x,c1=3.5,c2=6.75){
  zero <- -18.765
  time_scale <- 0.75
  out <- c1+c2*Bessel::AiryA(zero+qq(time_scale*abs(x)))
}

qq <- function(x){
  x - 0.017*x*x
}

#' Fit MCMC/Gibbs model
#' Mostly fits the Harvey and Trimbur model, except as noted here where some
#' parameters are more constrained
#' @param scale_A Scaling of amplitude
#' @param scale_Aj1 Scaling of amplitude response to Qr
#' @param scale_gamma Scaling of phase fluctuation
#' @export
jay_flinchem_chirptest <- function(scale_A=1.1,scale_Aj1=1.24,scale_gamma=1.){
  t = seq(-40.*24, 40.*24, 0.25)  # time in hours
  omega <- c(1.,2.,3.,4)
  omega <- c(1.,2.,3.,4)
  gamma <- scale_gamma*c(5.,7.5,15.,15.)*2.*pi/360.
  D1 <- tide_freq_radhr[['M1']]
  tnorm <- t*D1
  Qr <- Qr_reflect(t/24.)

  output <- list(time=t,Qr=Qr)

  A <- c(1.,2.,0.39,0.18)
  A <- c(0.75,1.45,0.33,0.15)*scale_A
  Aj0 <- A
  Aj1 <- c(0.33,0.33,0.25,0.2)*scale_Aj1
  Aj1 <- c(0.31,0.3,0.21,0.2)*scale_Aj1


  D <- Qr*1.

  for (i in seq(1,4)){
    name <- paste0("D",omega[i])
    phij <- gamma[i]*sqrt(Qr-1)
    Aj <- Aj0[i]*(1.-Aj1[i]*sqrt(Qr))
    constituent <- Aj*cos(tnorm*omega[i]-phij)
    output[[name]] <- constituent
    output[[paste0(name,"_amp") ]] <- Aj
    output[[paste0(name,"_phase") ]] <- phij
    output[[paste0(name,"_phasearg") ]] <- tnorm*omega[i]-phij
    D <- D + constituent
  }
  output[["D"]] <- D

  return(output)
}
