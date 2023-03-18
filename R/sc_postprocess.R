
modefunc <- function(x){
  tabresult <- tabulate(x)
  if(sum(tabresult == max(tabresult))>1) themode <- NA
  return(themode)
}

#' Take the mode of every column
#'
#' @param mat1 matrix over which modes are taken
#' @return vector of modes
#' @export
colModes <- function(mat1){
  apply(mat1, 2, modefunc)
}

#' Extract params from Gibbs sc_out output
#'
#' Constructs the matrices \code{FF} and \code{GG} used by DLM to represent
#' the stochastic cycle plus trend model
#' @param gibbs_output output of sc_gibbs
#' @param nburn number of cycles to discard as MCMC burn-in. If <0, first half is discarded
#' @param filter_eval if not null, should be a list including filter_freqs, filter_targets and filter_scale
#' @return sc_params Posterior mean parameters from model
#' @export
sc_params_from_gibbs <- function(gibbs_out,nburn=-1,filter_eval=NULL){
  g_sigma2_eps <- gibbs_out$sigma2_epsilon
  ntime <- length(g_sigma2_eps)
  if (nburn < 0){
    nburn = ntime %/% 2
  }
  n_iter <- length(g_sigma2_eps)
  sigma2_eps <- mean(g_sigma2_eps[(nburn+1):n_iter])
  sigma2_zeta <- mean(gibbs_out$sigma2_zeta[(nburn+1):n_iter])
  sigma2_kappa <- colMeans(gibbs_out$sigma2_kappa[(nburn+1):n_iter,])

  if ("rho" %in% names(gibbs_out)){
    rho <- mean(gibbs_out$rho[(nburn+1):n_iter])
  }else{
    rho <- 1.0 # unit root. #todo: make this easier
  }
  if(!is.null(filter_eval)){
    order_trend <- gibbs_out$order_trend
    order_cycle <- gibbs_out$order_cycle
    freqs <- gibbs_out$freqs

    filter_scale <- filter_eval$filter_scale
    filter_freqs <- filter_eval$filter_prior_freqs
    test_targets <- filter_eval$test_targets
    prior_val <- c()
    for (iter in (nburn+1):n_iter){
        eval_eps <- g_sigma2_eps[iter]
        eval_rho <- sc_out$rho[iter]
        qzeta <- sc_out$rho[iter]/eval_eps
        qkappa <- sc_out$sigma2_kappa[iter,]/eval_eps
        prior_val_new <- filter_prior(filter_freqs,test_targets,filter_scale,
                     gibbs_out$freqs,order_trend,order_cycle,rho,qzeta,qkappa)
        prior_val <- c(prior_val,prior_val_new)

    }
    prior_val <- mean(prior_val)


  }else{
    prior_val <- NA
  }

  freqs <- gibbs_out$freqs
  order_trend <- gibbs_out$order_trend
  order_cycle <- gibbs_out$order_cycle
  sc_params <- list(sigma2_eps=sigma2_eps,sigma2_zeta=sigma2_zeta,sigma2_kappa=sigma2_kappa,
                    rho=rho,freqs=freqs,order_trend=order_trend,order_cycle=order_cycle,
                    filter_prior_val=prior_val)
}



#' Build the stochastic cycle model plus trend matrices from Gibbs sc_out
#'
#' Constructs the matrices \code{FF} and \code{GG} used by DLM to represent
#' the stochastic cycle plus trend model
#' @param gibbs_output output of sc_gibbs
#' @param nburn number of cycles to discard as MCMC burn-in. If <0, first half is discarded
#' @return model dlm model
#' @export
sc_model_from_gibbs <- function(gibbs_out,nburn=-1,sigma2_diffuse=1e3){
  parms <- sc_params_from_gibbs(gibbs_out,nburn=nburn)
  print(parms)
  model <- sc_model_build(parms$order_trend,parms$order_cycle,
                          parms$freqs,parms$rho,
                          parms$sigma2_eps,parms$sigma2_zeta,
                          parms$sigma2_kappa,sigma2_diffuse)
  model
}


#' Constructs the matrices \code{FF} and \code{GG} used by DLM to represent
#' the stochastic cycle plus trend model
#' @param parms list of parameters, possibly output of sc_params_from_gibbs
#' @return model dlm model
#' @export
sc_model_from_parms <- function(parms){
  model <- sc_model_build(parms$order_trend,parms$order_cycle,parms$freqs,parms$rho,
                          parms$sigma2_eps,parms$sigma2_zeta,parms$sigma2_kappa)
  model
}




#' Build the stochastic cycle model plus trend matrices
#'
#' Constructs the matrices \code{FF} and \code{GG} used by DLM to represent
#' the stochastic cycle plus trend model
#' @param gibbs_output output of sc_gibbs
#' @param nburn number of cycles to discard as MCMC burn-in. If <0, first half is discarded
#' @return model dlm model
#' @export
sc_model_from_gibbs_mode <- function(gibbs_out,nburn=-1){
  g_sigma2_eps <- gibbs_out$sigma2_epsilon
  ntime <- length(g_sigma2_eps)
  if (nburn < 0){
    nburn = ntime %/% 2
  }
  n_iter <- length(g_sigma2_eps)
  sigma2_eps <- mode(g_sigma2_eps[(nburn+1):n_iter])
  sigma2_zeta <- mode(gibbs_out$sigma2_zeta[(nburn+1):n_iter])
  sigma2_kappa <- colModes(gibbs_out$sigma2_kappa[(nburn+1):n_iter,])
  print("sigma_kappa")
  print(sigma2_kappa)
  if ("rho" %in% names(gibbs_out)){
      rho <- mode(gibbs_out$rho[(nburn+1):n_iter])
  }else{
    rho <- 1. # unit root
  }
  freqs <- gibbs_out$freqs
  order_trend <- gibbs_out$order_trend
  order_cycle <- gibbs_out$order_cycle

  model <- sc_model_build(order_trend,order_cycle,freqs,rho,sigma2_eps,sigma2_zeta,sigma2_kappa)
  model
}


#' recreate full state
#' @param state s
#' @param order_trend order of trend
#' @param order_cycle order of cycle
#' @param nfreq number of frequencies
#' @export
reconstruct_state <- function(state,order_trend,order_cycle,nfreq)
{
  subtide <- state[,1]
  diurnal <- state[,1+order_trend]
  recon <- subtide
  for (icycle in 1:nfreq){
    recon <- recon + state[,order_trend + (icycle-1)*2*order_cycle + 1]
  }
  list(subtide=subtide,full=recon,diurnal=diurnal)

}

#' do plot
#' @param ysave saved data in case y has gaps
#' @param y data with gaps
#' @param recon reconstructed state
#' @param select selection to plot
#' @export
do_plot <- function(ysave,y,recon,select,harmonic=0.){
  par(mfrow=c(2,1),mar=c(2,2,2,2))

  ts.plot(ysave[select],col="gray",lwd=3)
  lines(y[select],col="red")
  lines(recon$full[select],col="darkgreen")
  lines(recon$subtide[select],col="blue")
  if (is.vector(harmonic)){
    full = recon$full[select] + harmonic[select]
  }else{
    full = recon$full[select]
  }
  lines(full,col="dark green")
  #ts.plot(recon$diurnal)
}

#' Extract amplitudes, phase and values from state output matrix
#' @param state state to interpret
#' @param order_trend order of trend
#' @param order_cycle order of cycle
#' @param nfreq number of frequencies
#' @param flabels names for freqs
#' @export
interpret_state <- function(state,order_trend,order_cycle,nfreq,flabels){
  fitted  <- reconstruct_state(state,order_trend,order_cycle,nfreq)
  subtide <-  state[,1]
  freq <- list()
  amplitude <- list()
  phase <- list()
  for (jj in (1:nfreq)) {
    # todo: hardwired
    label <- flabels[jj]
    jndx = order_trend + (jj-1)*2*order_cycle + 1
    amplitude[[label]] <- sqrt(rowSums(state[,jndx:(jndx+1)]^2.))
    phase[[label]] <- atan2(state[,jndx+1],state[,jndx])
    freq[[label]] <- state[,jndx]
  }
  list(fitted=fitted,subtide=subtide,freq=freq,amplitude=amplitude,phase=phase)
}

quart <- function(x){
  quantile(x,probs=c(0.25,0.5,0.75))
}


#' Compute posterior mean after burnin
#' @param x data
#' @param nburn number of burnin to skip, less than zero means half
#' @param is_angle are the values angles (uses vector averaging)
#' @param quantiles vector of quantiles to extact
#' @export
aggregate_iter <- function(x,nburn=-1,is_angle=FALSE,quantiles=c()){
  if (is_angle & !(is.null(quantiles))){
    stop("quantiles not available for angles")
  }
  niter <- dim(x)[1]
  if (nburn < 0){
    nburn = niter %/% 2
  }
  print(nburn)
  print(niter)
  if (is_angle){
    xave <- colMeans(exp(1i*x[(nburn+1):niter,]))
    xave <- Arg(xave)
  }else{
    if (is.null(quantiles)){
      xave <- colMeans(x[(nburn+1):niter,])
    }else{
      xave = apply(x[(nburn+1):niter,],2,quart)
    }
  }
  return(xave)
}

#' Summarize the state info from Gibbs run and interpretation of run
#' @param sc_out Result of Bayes computation
#' @param interp Single pass result from interpret_state
#' @param out_dt Time step to use for labeling output (e.g. 0.25 for 15min)
#' @param nburn number of burnin to skip, less than zero means half
#' @param outfile Name of output file for resulting dataframe (or NULL)
#' @export
summarize_scout_result <- function(sc_out,interp,out_dt=0.25,nburn=10,outfile=NULL){

  flabels <- names(sc_out$amplitude_thinned)

  # This calculate the number of outputs, compensating for the extra value at t=0 that can
  # cause mismatch with the original series
  nout=length(interp$subtide)-1   # nout should be same length as original series
  nplusone <- length(interp$subtide)
  times <- seq(0.,by=out_dt,length.out=nout)

  # aggregate
  ip1a <- aggregate_iter(sc_out$subtide_thinned,nburn=nburn)
  ip1b <- aggregate_iter(sc_out$subtide_thinned,quantiles=c(0.,0.5,0.75),nburn=nburn)

  # Create a list that later will be converted to the output dataframe
  # starting with
  # parm means the state that results from the additional run wiht posterior mean parameters (interp)
  # postmean means the posterior mean of the state from MCMC process
  # 25, 50, 75 are the respective quantiles of the result
  archive = list(times=times,
                 subtide_parm=interp$subtide[2:nplusone],
                 subtide_postmean=ip1a,
                 subtide25= ip1b[1,],
                 subtide50 =ip1b[2,],
                 subtide75= ip1b[3,])

  # Loop through frequencies and add analogous quantities
  # for the values, amplitudes and phase (fewer for phase
  # since quantiles are a bit harder to uniquely define defined)
  nfreq = length(sc_freq)
  for (ifreq in 1:nfreq){
    flab0 = flabels[ifreq]
    flab = flabels[ifreq]

    print("Processing freq")
    # Value of frequency

    postmean <- aggregate_iter(sc_out$freq_thinned[[flab0]],nburn=nburn)
    anal <- aggregate_iter(sc_out$freq_thinned[[flab0]], nburn=nburn,quantiles=c(0.25,0.5,0.75))


    archive[[paste0(flab,'_val_parm')]] <- interp$freq[[flab0]][2:nplusone]
    archive[[paste0(flab,'_val_postmean')]] <- postmean[2:nplusone]
    archive[[paste0(flab,'_val_25')]] <- anal[1,]
    archive[[paste0(flab,'_val_50')]] <- anal[2,]
    archive[[paste0(flab,'_val_75')]] <- anal[3,]

    print("Processing amplitude")

    # Amplitude statistics
    postmean <- aggregate_iter(sc_out$amplitude_thinned[[flab0]],nburn=nburn)
    anal <- aggregate_iter(sc_out$amplitude_thinned[[flab0]],nburn=nburn,quantiles=c(0.25,0.5,0.75))
    archive[[paste0(flab,'_amp_parm')]] <- interp$amplitude[[flab0]][2:nplusone]
    archive[[paste0(flab,'_amp_postmean')]] <- postmean[2:nplusone]
    archive[[paste0(flab,'_amp_25')]] <- anal[1,]
    archive[[paste0(flab,'_amp_50')]] <- anal[2,]
    archive[[paste0(flab,'_amp_75')]] <- anal[3,]

    # Phase statistics
    postmean <- aggregate_iter(sc_out$phase_thinned[[flab0]],nburn=nburn,is_angle=TRUE)
    archive[[paste0(flab,'_phase_parm')]] <- interp$phase[[flab0]][2:nplusone]    # ws ip1a
    archive[[paste0(flab,'_phase_postmean')]] <- postmean[2:nplusone]
  }

  archive_df <- as.data.frame(do.call(cbind,archive))
  row.names(archive_df) <- format(archive_df[,'times'],digits=2)

  if (!is.null(outfile)){
    write.csv(format(archive_df,digits=5),file=outfile,quote=FALSE,row.names=FALSE)
  }

  archive_df
}


