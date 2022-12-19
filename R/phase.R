
#' Remove linear phase
#' @param result phase from which to remove reference
#' @param freq frequency for time-dependent part of phase
#' @param times times for time-dependent Part of phase
#' @return phase relative to linear phase
#' @export
rm_linear_phase<-function(result,freq,times){
  if(is.null(dim(result))){
    uresult <- signal::unwrap(result)
    nt <- length(result)
    nsim = 1
    inphase <- signal::unwrap(times*freq)

  }else{
    uresult <- t(apply(result,1,signal::unwrap))
    nt <-  dim(uresult)[2]
    nsim <- dim(uresult)[1]
    inphase <- times*freq
    inphase <- matrix(data=times*freq,nrow=nsim,ncol=nt,byrow=TRUE)
  }
  rmphase <- uresult - inphase
  return(rmphase)
}

# 'Find multiple of 2$\pi$ closest to the median of the input vector x
# ' @param x vector of values
med_mod2pi <- function(x){
  med <- median(x)
  twopi <- 2.*pi
  round(med/(twopi),digits=0)*twopi
}

#' Removes a reference phase argument from a phase result
#'
#' Comparing phases is tricky because of the wrapping and
#' other modulo \eqn{2\pi} problems, but also because in some test cases the frequency
#' may be inexact which causes issues over long spans of time.
#' This function compares the total phase argument including the linear in time
#' part (\eqn{\omega t + \phi}). After unwrapping/comparing the modulo \eqn{2\pi}
#' issue is addressed by calculating the median residual and
#' removing multiple of \eqn{2\pi} that is closest.
#' @param result result phase from which to remove reference
#' @param freqarg reference argument including freq*t part and perturbation
#' @param phi additive scalar or time-dimension reference to be added
#' @return phase relative to linear phase
#' @export
rm_reference_phase<-function(result,freqarg){
  if(is.null(dim(result))){
    print("vec")
    uresult <- signal::unwrap(result)
    fresult <- signal::unwrap(freqarg)
    rmphase <- uresult - fresult

    closest_mult2pi = med_mod2pi(rmphase)
    rmphase <- rmphase - closest_mult2pi
    return(rmphase)

  }else{
    print("mat")
    print(dim(result))
    uresult <- apply(result,1,signal::unwrap) # transposes so time is rows
    fresult <- signal::unwrap(freqarg)
    rmphase <- uresult - fresult
    rmphase <- sweep(rmphase,2,apply(rmphase,2,med_mod2pi))
    return(t(rmphase))
  }
}
