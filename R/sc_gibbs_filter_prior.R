

#' Fit MCMC/Gibbs model using the Harvey and Trimbur model with known frequency
#' @param y The observed data
#' @param order_trend order of the trend (low pass) component of the model
#' @param order_cycle order of the cyclical part of the model
#' @param freqs frequencies used in fitting
#' @param samples_per_hr number of obs per hour
#' @param mc_iter number of iterations
#' @param filter_prior_freqs filter quality test freqs in cpd
#' @param test_targets filter quality test targets (gain)
#' @param filter_scale scale of loss/prior on bad filter fit to freqs/targets
#' @param initial_sigma2_epsilon initial MCMC iterate for sigma2_epsilon
#' @param mh_sd_eps Metropolis Hastings random walk std. deviation for sigma2_epsilon
#' @param sigma2_zeta_fixed is sigma_zeta_fixed (and thus out of MCMC)?
#' @param sigma2_diffuse magnitude for diffuse initial condition
#' @param initial_sigma2_zeta initial MCMC value for sigma2_zeta
#' @param mh_sd_zeta Metropolis Hastings random walk std. deviation for zeta
#' @param initial_sigma2_kappa initial variance for kappa, either scalar or one per freq
#' @param mh_sd_kappa Metropolis Hastings random walk std. deviation for kappa, either scalar or one per freq
#' @param rho_fixed is rho fixed (and thus out of MCMC)
#' @param initial_rho starting rho
#' @param rho_high top of the uniform prior, often guided by conditioning of problem
#' @param rho_low bottom of the uniform prior, often guided by need for neat partiion
#' @param mh_sd_rho Metropolis Hastings random walk std. deviation for rho
#' @param thin_rate thinning rate for posterior draws (every thin_rate sample is used)
#' @export
sc_gibbs_filter_prior <- function(y,order_trend,order_cycle,freqs,samples_per_hr,
                                  mc_iter,
                                  filter_prior_freqs,test_targets,filter_scale,
                                  initial_sigma2_epsilon=1.e-4,mh_sd_eps=NULL,
                                  sigma2_zeta_fixed=FALSE,sigma2_diffuse=1.e4,
                                  initial_sigma2_zeta=1.e-6,mh_sd_zeta=NULL,
                                  initial_sigma2_kappa=1.e-11,mh_sd_kappa=NULL,
                                  rho_fixed=FALSE,initial_rho=0.9925,
                                  rho_low=0.985,
                                  rho_high=0.9925,
                                  mh_sd_rho=NULL,
                                  thin_rate=5)
{


  # Preliminaries constructing model and doing regression if it is offline
  nobs <- length(y) # total length of obs data, as T in the equation
  miss <- is.na(y)
  ygood <- y[!miss]
  ngood <- length(ygood)

  # Set Metropolis Hastings jump std. deviation if not provided
  if (is.null(mh_sd_eps)){
    mh_sd_eps = initial_sigma2_epsilon/15.
  }
  if (is.null(mh_sd_zeta)){
    mh_sd_zeta = initial_sigma2_zeta/15.
  }
  if (is.null(mh_sd_kappa)){
    mh_sd_kappa = initial_sigma2_kappa/15.
  }
  if (is.null(mh_sd_rho)){
    mh_sd_rho = initial_rho/15.
  }

  n_freq <- length(freqs)
  # todo: this seems not parallel to the sc_gibbs function
  sc_freq <- freqs*samples_per_hr
  rho <-  initial_rho
  scmat <- sc_model_mat(order_trend, order_cycle, freqs, rho)

  V <- diag(c(1))

  # W is state innovation variance,
  Wtrend <- c(rep(0.,order_trend-1),1.)*initial_sigma2_zeta
  if (length(initial_sigma2_kappa) == 1){
    initial_sigma2_kappa <- rep(initial_sigma2_kappa,n_freq)
  }
  # todo: Moved this out of brackets above because it seemed wrong
  # does this work for length > 1 and ==1 ?
  Wcycle <- rep( c(rep(0,2*(order_cycle-1)),1,1),n_freq)*initial_sigma2_kappa

  if(length(mh_sd_kappa)==1){
    mh_sd_kappa = rep(mh_sd_kappa,n_freq)
  }

  mat_zero <- matrix(0, ncol = n_freq, nrow = 2*(order_cycle-1))
  Wcycle <- as.vector(matrix( rbind(mat_zero,initial_sigma2_kappa,initial_sigma2_kappa),byrow=TRUE))

  W <- diag(c(Wtrend,Wcycle))


  nstate <- order_trend + 2*order_cycle*n_freq

  # calculate covariance matrix for nth-order cycle
  #covar_matrix_size = order_cycle**2
  #These are initial calculations that are updated as rho evolves
  covar_matrix_cycle <- list()
  covar_matrix_cycle_inv <- list()
  covar_matrix_multiplied <- list()

  # These are to cache some covariance calculations for proposed new rho that might
  # not be accepted
  covar_new <- list()
  covar_inv_new <- list()


  for (jj in (1:n_freq)) {
    covar_matrix_cycle[[jj]] <- calc_cycle_covariance_matrix(order_cycle,rho,freqs[jj])
    covar_matrix_cycle_inv[[jj]] <- solve(covar_matrix_cycle[[jj]])
    covar_matrix_multiplied[[jj]] <- covar_matrix_cycle[[jj]]*initial_sigma2_kappa[[jj]]
  }

  # Initial condition for the state is zero but evolves (check this)
  m0 <- rep(0, dim(W)[1])
  m0[1] <- 0.

  # Uninformative C0
  C0_trend <- diag(1,order_trend)*sigma2_diffuse

  C0_cycle <- dlm::bdiag(covar_matrix_multiplied)
  C0 <- dlm::bdiag(C0_trend, C0_cycle)
  model <- dlm::dlm(FF=scmat$FF,V=V,GG=scmat$GG,W=W,m0=m0,C0=C0)

  # begin MCMC/Gibbs' sampling iteration
  print (paste('Start MCMC',date()))
  print (paste('Requested # MCMC iterations (mc_iter) =',mc_iter))
  flush.console()

  # These are for logging samples obtained during MCMC process
  #beta_draw_save <- matrix(numeric(mc_iter*k), nrow = mc_iter, ncol = k) # todo: removed intercept so k not k+1
  sigma2_epsilon_save <- numeric(mc_iter)
  sigma2_kappa_save <- matrix(numeric(mc_iter*n_freq), nrow = mc_iter, ncol = n_freq)
  sigma2_zeta_save <- numeric(mc_iter)
  rho_save <- matrix(numeric(mc_iter),nrow=mc_iter,ncol=1)
  accept_reject <- list(eps=c(0,0),zeta=c(0,0),
                        kappa=lapply(freqs,function(x){c(0,0)}),rho=c(0,0))
  subtide_thinned <- c()
  amplitudes_thinned <- list()
  phase_thinned <- list()
  freq_thinned <- list()
  state_thinned <- c()

  alltargs <- regular_filter_targets()
  testfreqs <- alltargs$filter_prior_freqs
  test_targets <- alltargs$test_targets


  sigma2_zeta <- initial_sigma2_zeta
  sigma2_kappa <- rep(initial_sigma2_kappa,n_freq)
  sigma2_eps_draw <- initial_sigma2_epsilon




  set.seed(201)
  for (itr in(1:mc_iter)) {
    print("========")
    print(paste('MCMC iteration',itr))

    #step 1 - sample beta (normal) and sigma2 (inverse gamma)
    if (itr==1) {
      # todo, eliminated the mean subtraction
      y_mod <- y # - mean(y,na.rm=TRUE)
    } else {
      state_part <-alpha_draw %*% t(model$FF)
      y_mod <- y - state_part
    }

    tide_res = y_mod[!miss]    # not really a tide residual in this case

    s2 <- crossprod(tide_res,tide_res) #equation (14.7) in BDA book
    dim(s2) <- NULL


    eps_old <- sigma2_eps_draw
    eps_new <- eps_old + rnorm(1,mean=0,sd=mh_sd_eps)
    dense_eps_old <- -0.5*nobs*log(eps_old)-0.5*s2/eps_old
    if (eps_new>0.){
      dense_eps_new <- -0.5*nobs*log(eps_new)-0.5*s2/eps_new
    }else{
      dense_eps_new = 0.
    }
    qzeta <- sigma2_zeta/eps_old
    qkappa <- sigma2_kappa/eps_old
    cost_old <- filter_prior(filter_prior_freqs=testfreqs,test_targets=test_targets,scale=filter_scale,
                             sc_freq=sc_freq,m=order_trend,n=order_cycle,
                             rho=rho,qzeta=qzeta,qkappa=qkappa)
    qzeta <- sigma2_zeta/eps_new
    qkappa <- sigma2_kappa/eps_new
    cost_new <-filter_prior(testfreqs,test_targets,scale=filter_scale,
                            sc_freq,m=order_trend,n=order_cycle,
                            rho=rho,qzeta=qzeta,qkappa=qkappa)
    accept_ratio_eps = exp(-cost_new + cost_old + dense_eps_new - dense_eps_old)
    u <- runif(1,0.,1.)
    accept_eps <- (u <= accept_ratio_eps)
    if (accept_eps){
      sigma2_eps_draw <- eps_new
      accept_reject$eps[1] <- accept_reject$eps[1] + 1
    }else{
      accept_reject$eps[2] <- accept_reject$eps[2] + 1
    }
    print("Epsilon")
    print(paste("sigma2_eps_old",eps_old,"sigma2_eps_new",eps_new))
    print(paste("prior-based cost_new",cost_new,"cost_old",cost_old))
    print(paste("sigma2_eps accepted: ",accept_eps))
    print(paste("Current sigma2_eps=",sigma2_eps_draw," scale version:",s2/(nobs)))

    dlm::V(model) <- matrix(sigma2_eps_draw)
    sigma2_epsilon_save[itr] <- sigma2_eps_draw


    # Have to embed these back in the full size including missing data
    res_good <- ygood
    res <- y # for size including missing data, don't care about y
    res [miss] <- NA
    res[!miss] <- res_good

    #step 2 - sample state vector alpha using dlmBSample
    modFilt <- dlm::dlmFilter(res, model)

    alpha_draw <- dlm::dlmBSample(modFilt)
    alpha_draw <- alpha_draw[-1,]
    m0_new <- alpha_draw[1,]
    state_part2 <-alpha_draw %*% t(model$FF)
    do_thinning <- (itr %% thin_rate) == 0
    if (do_thinning){
      reconstruct <- state_part2
      print("finished state draw")
    }

    #step 3 - sample sigma_zeta (trend) and sigma_kappa (cycle),
    #         both follow inverse gamma distribution
    # notation in Harvey and Trimbur is beta(n) - beta(n-1) but we already used beta
    trend_diff <- alpha_draw[(2:nobs),order_trend]-alpha_draw[1:(nobs-1),order_trend]
    scale_zeta <- crossprod(trend_diff,trend_diff)
    shape_zeta <- nobs/2+.5

    print(paste("shape_zeta=",shape_zeta,"scale_zeta=",scale_zeta))
    if (sigma2_zeta_fixed){
      sigma2_zeta <- initial_sigma2_zeta
      sigma2_zeta_fixed <- FALSE
    }else{
      zeta_old <- sigma2_zeta
      zeta_new <- zeta_old + rnorm(1,mean=0,sd=mh_sd_zeta)
      dense_zeta_old <- -0.5*(nobs-1)*log(zeta_old)-0.5*scale_zeta/zeta_old
      if (zeta_new>0.){
       dense_zeta_new <- -0.5*(nobs-1)*log(zeta_new)-0.5*scale_zeta/zeta_new
      }else{
        dense_zeta_new <- 0.
      }
      qzeta <- zeta_old/sigma2_eps_draw
      qkappa <- sigma2_kappa/sigma2_eps_draw
      cost_old <- filter_prior(filter_prior_freqs=testfreqs,test_targets=test_targets,scale=filter_scale,
                          sc_freq=sc_freq,m=order_trend,n=order_cycle,
                          rho=rho,qzeta=qzeta,qkappa=qkappa)
      qzeta <- zeta_new/sigma2_eps_draw
      cost_new <-filter_prior(testfreqs,test_targets,scale=filter_scale,
                              sc_freq,order_trend,order_cycle,rho,qzeta,qkappa)
      accept_ratio_zeta = exp(-cost_new + cost_old + dense_zeta_new - dense_zeta_old)
      u <- runif(1,0.,1.)
      accept_zeta <- (u <= accept_ratio_zeta)
      if (accept_zeta){
        sigma2_zeta <- zeta_new
        accept_reject$zeta[1] <- accept_reject$zeta[1] + 1
      }else{
        accept_reject$zeta[2] <- accept_reject$zeta[2] + 1
      }
    }
    print("Zeta")
    print(paste("sigma2_zeta_old",zeta_old,"sigma2_zeta_new",zeta_new))
    print(paste("cost_new",cost_new,"cost_old",cost_old))
    print(paste("acceptance rate",accept_ratio_zeta,"draw",u))
    print(paste("sigma2_zeta accepted: ",accept_zeta))
    print(paste("sigma2_zeta=",sigma2_zeta," scale version:",scale_zeta/(nobs)))
    sigma2_zeta_save[itr] <- sigma2_zeta
    Wtrend <- c(rep(0.,order_trend-1),sigma2_zeta)
    print("*-*")

    #3.2 - work on cycle part
    Wcycle <- c()
    C0_cycle_new <- list()
    rho_sel <- rho  # todo: Too many copies ... relic of them being independent
    rho_old <- rho_sel
    # There is only one value of rho for all freqs
    rho_new <- rho_old + rnorm(1,mean=0,sd=mh_sd_rho)
    # These are simple prob bounds, not a prior like in other versions
    rho_in_bounds <- ((rho_new >= rho_low) & (rho_new <= rho_high))

    # Accumulates log probability for rho, which requires summing across all
    # Frequencies
    pold <- 0
    pnew <- 0

    print("Kappa")
    for (jj in (1:n_freq)) {
      # each freq has its own variance
      freqs_sel <- freqs[jj]
      freq_name <- names(freqs)[jj]
      covar_matrix_inv_sel <- covar_matrix_cycle_inv[[jj]]
      scale_kappa <- cal_shape_kappa(order_trend, order_cycle, jj, freqs_sel,
                                     rho_sel, covar_matrix_inv_sel,
                                     nobs, alpha_draw)
      # todo: review scale
      shape_kappa <- nobs+order_cycle-1
      #sigma2_kappa <- invgamma::rinvgamma(1, rate=scale_kappa/2., shape=shape_kappa)

      kappa_old <- sigma2_kappa
      kappa_new <- sigma2_kappa
      kappa_new[[jj]] <- kappa_old[[jj]] + rnorm(1,mean=0,sd=mh_sd_kappa[jj])

      dense_kappa_old <- -(nobs+order_cycle-1)*log(kappa_old[[jj]])-0.5*scale_kappa/kappa_old[[jj]]
      if(kappa_new[[jj]]>0){
        dense_kappa_new <- -(nobs+order_cycle-1)*log(kappa_new[[jj]])-0.5*scale_kappa/kappa_new[[jj]]
      }else{
        dense_kappa_new = 0.0
      }
      qzeta <- sigma2_zeta/sigma2_eps_draw
      qkappa <- kappa_old/sigma2_eps_draw
      cost_old <-filter_prior(testfreqs,test_targets,scale=filter_scale,
                              sc_freq,order_trend,order_cycle,rho,qzeta,qkappa)
      qkappa <- kappa_new/sigma2_eps_draw
      cost_new <-filter_prior(testfreqs,test_targets,scale=filter_scale,
                              sc_freq,order_trend,order_cycle,rho,qzeta,qkappa)
      print(paste("filter prior cost_old",cost_old,"cost_new",cost_new,"dense_kappa_old",dense_kappa_old,
                  "dense_kappa_new",dense_kappa_new))

      accept_ratio_kappa = exp(-cost_new + dense_kappa_new +cost_old - dense_kappa_old)
      u <- runif(1,0.,1.)
      accept_kappa <- (u <= accept_ratio_kappa)
      print(paste("accept ratio",accept_ratio_kappa))
      if (accept_kappa){
        sigma2_kappa[[jj]] <- kappa_new[[jj]]
        accept_reject$kappa[[freq_name]][1] <- accept_reject$kappa[[freq_name]][1] + 1
      }else{
        accept_reject$kappa[[freq_name]][2] <- accept_reject$kappa[[freq_name]][2] + 1
      }

      print(paste("accepted sigma2 in index",jj,accept_kappa))
      print(paste("freq jj =",jj," shape=",shape_kappa," scale=",scale_kappa,
                  " sigma2=",sigma2_kappa[[jj]]))

      sigma2_kappa_save[itr,jj] <- sigma2_kappa[[jj]]
      Wcycle <- c(Wcycle,rep(0,2*(order_cycle-1)),sigma2_kappa[[jj]],sigma2_kappa[[jj]])

      if ((!rho_fixed) & rho_in_bounds){
        # covar_matrix_cycle only depends on rho and freq, so this is still
        # ok, unperturbed by work on zeta and eps
        covar_old<- covar_matrix_cycle[[jj]]
        covar_inv_old <- covar_matrix_cycle_inv[[jj]]

        arg_old <- scale_kappa
        arg_old <- arg_old/sigma2_kappa[[jj]]
        pdivold <- 0.5*(determinant(sigma2_kappa[[jj]]*covar_old,
                                    logarithm=TRUE)$modulus[[1]])

        pold <- pold -0.5*arg_old-pdivold   # this is the log  probability

        # now rho is changed, so recalculate
        covar_new[[jj]] <- calc_cycle_covariance_matrix(order_cycle,rho_new,freqs[jj])
        solve_out <- tryCatch(solve(covar_new[[jj]]),error=function(c){NULL})
        if (is.null(solve_out)){
          print("Discarding new value of rho due to covariance conditioning")
          # This is meant to be so large and negative that it doesn't get used
          pnew <- -1000000.
          rho_in_bounds = FALSE
        }else{
          covar_inv_new[[jj]] <- solve_out
          arg_new <- cal_shape_kappa(order_trend, order_cycle, jj, freqs_sel,
                                   rho_new, covar_inv_new[[jj]], nobs, alpha_draw)
          arg_new <- arg_new/sigma2_kappa[[jj]]
          pdivnew <- 0.5*(determinant(sigma2_kappa[[jj]]*covar_new[[jj]],logarithm=TRUE)$modulus[[1]])
          pnew <- pnew -0.5*arg_new-pdivnew   # this is the log probability
        }
      }
    }  # loop jj over frequencies

    # Now calculate prior, which is determined across all freqs at once
    qzeta <- sigma2_zeta/sigma2_eps_draw
    qkappa <- sigma2_kappa/sigma2_eps_draw

    # Already on log scale. Must be multiplied by -1 to go from cost to prob
    prior_old <-filter_prior(testfreqs,test_targets,scale=filter_scale,
                            sc_freq,order_trend,order_cycle,rho_old,qzeta,qkappa)
    prior_new <-filter_prior(testfreqs,test_targets,scale=filter_scale,
                            sc_freq,order_trend,order_cycle,rho_new,qzeta,qkappa)

    print(paste("Rho prior_old",prior_old,"prior_new",prior_new))
    pold <- pold - prior_old
    pnew <- pnew - prior_new

    # Now do the M-H step for rho
    # Initialize with rho_old, assuming no step overwrite with rho_new if accepted
    rho <- rho_old    # This doesn't do anything, but it clarifies
    # If the new value is precluded by being out of bounds no action needed except to record it
    if ((!rho_fixed) & rho_in_bounds){
      if (pnew >= pold){
        # ratio > 1, accept without transformations and random draw
        print(paste("rho: accept ratio >= 1, new rho = ",rho_new))
        rho <- rho_new

      }else{
        logpdiff = pnew-pold
        accept_rate <- exp(logpdiff)
        u =runif(1,0.,1.)
        print(paste("rho: calc accept_rate=",accept_rate,"draw",u))
        if (u < accept_rate){
          print(paste("accept rho = ",rho_new))
          rho <- rho_new
        }else{
          print(paste("rho",rho_new,"rejected based on draw, retain ",rho_old))
        }
      } # pnew >= pold
    }else{
      # rho: rejected step because the new row went outside bounds or because it is fixed
      if (rho_fixed){
        rho <- rho_old
        print(paste("rho fixed, retain",rho_old))
      }else{
        rho <-rho_old
        print(paste("rho: reject b/c draw out of bounds (rho_low,rho_high)",
                     rho_new,", retain ",rho_old))
      }
    }

    if (rho != rho_old){
      if (rho !=  rho_new) { print("BUG in rho code")}
      rho_save[itr] <- rho_new
      accept_reject$rho[1] <- accept_reject$rho[1] + 1
      for (jj in 1:n_freq){
        covar_matrix_cycle[[jj]] <- covar_new[[jj]]
        covar_matrix_cycle_inv[[jj]] <- covar_inv_new[[jj]]
        C0_cycle_new[[jj]] <- covar_matrix_cycle[[jj]]*sigma2_kappa_save[[itr,jj]]
      }
    }else{
      rho <- rho_old # for clarity
      rho_save[itr] <- rho_old
      accept_reject$rho[2] <- accept_reject$rho[2] + 1
      for (jj in 1:n_freq){
        # todo: This is probably avoidable, but need to make sure C0_cycle_new is set for
        # the first iterations before a rho is accepted
        C0_cycle_new[[jj]] <- covar_matrix_cycle[[jj]]*sigma2_kappa_save[[itr,jj]]
      }
    }

    print("Accept/reject rates:")
    print(paste("epsilon :",paste(accept_reject$eps,collapse=" ")))
    print(paste("zeta: ",paste(accept_reject$zeta,collapse=" ")))
    print(paste("rho",paste(accept_reject$rho,collapse=" ")))
    for (jj in 1:n_freq){
      freq_name <- names(freqs)[jj]
      arrate <- accept_reject$kappa[[freq_name]]
      print(paste("kappa[",freq_name,"]",paste(arrate,collapse=" ")))
    }
    print("*")


    dlm::W(model) <- diag(c(Wtrend,Wcycle))
    C0_cycle <- dlm::bdiag(C0_cycle_new)
    dlm::C0(model) <- dlm::bdiag(C0_trend, C0_cycle)

    if (do_thinning){
      print(paste("iter=",itr))
      print ("thinned")
      state_thinned <- rbind(state_thinned,t(reconstruct))
      subtide_thinned <-  rbind(subtide_thinned,t(alpha_draw[,1]))
      for (jj in (1:n_freq)) {
        flabels <- names(freqs)
        label <- flabels[jj]
        label0 <- paste(flabels[jj],"_cos",sep="")
        label1 <- paste(flabels[jj],"_sin",sep="")
        jndx = order_trend + (jj-1)*2*order_cycle + 1
        amplitudes_thinned[[label]] <- rbind(amplitudes_thinned[[label]],t(sqrt(rowSums(alpha_draw[,jndx:(jndx+1)]^2.))))
        phase_thinned[[label]] <- rbind(phase_thinned[[label]],t(atan2(alpha_draw[,jndx+1],alpha_draw[,jndx])))
        freq_thinned[[label]] <- rbind(freq_thinned[[label]],t(alpha_draw[,jndx]))
      }
    }
  }

  print (paste('Finish MCMC',date()))
  flush.console()
  return(list(beta=beta_draw_save,sigma2_epsilon=sigma2_epsilon_save,sigma2_zeta=sigma2_zeta_save,
              sigma2_kappa = sigma2_kappa_save,rho = rho_save,
              order_trend = order_trend, order_cycle = order_cycle, freqs = freqs,
              state_thinned=state_thinned,
              subtide_thinned = subtide_thinned,
              amplitude_thinned=amplitudes_thinned,
              phase_thinned=phase_thinned,freq_thinned=freq_thinned,accept_reject=accept_reject))
}
