
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


#' Removes a reference phase argument including freq*t and perturbation phi(t)
#' @param result result phase from which to remove reference
#' @param phi perturbation part of reference
#' @param freqarg reference argument including freq*t part and perturbation
#' @return phase relative to linear phase
#' @export
rm_reference_phase<-function(result,phi,freqarg){
  if(is.null(dim(result))){
    print("vec")
    uresult <- signal::unwrap(result)
    fresult <- signal::unwrap(freqarg)
    rmphase <- uresult - fresult
    rmphase <- rmphase - median(rmphase) + phi
    return(rmphase)

  }else{
    print(dim(result))
    uresult <- apply(result,1,signal::unwrap) # transposes so time is rows
    fresult <- signal::unwrap(freqarg)
    rmphase <- uresult - fresult
    rmphase <- sweep(rmphase,2,apply(rmphase,2,median)) + phi
    return(t(rmphase))
  }
}
