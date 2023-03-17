#' Fit MCMC/Gibbs model
#' Mostly fits the Harvey and Trimbur model, except as noted here where some
#' parameters are more constrained
#' @param y The observed data
#' @param order_trend order of the trend (low pass) component of the model
#' @param order_cycle order of the cyclical part of the model
#' @param freqs frequencies used in fitting
#' @param mc_iter number of MCMC iterations to perform
#' @param initial_sigma2_epsilon initial condition for MCMC for sigma2_epsilon
#' @param sigma2_zeta_fixed treat sigma2_zeta as fixed
#' @param initial_sigma2_zeta initial condition for MCMC for sigma2_zeta. See defaults.mcmc
#' @param initial_sigma2_kappa initial condtion for MCMC for sigma2_kappa
#' @param rho_fixed whether to fix rho rather than include in MCMC
#' @param initial_rho initial rho for MCMC, also used for rho if rho_fixed=TRUE
#' @param rho_high upper bound on rho, if same as rho_low to 5 digits rho will be fixed
#' @param rho_low lower bound on rho uniform. Initial rho is average of rho_high and rho_low
#' @param mh_sd_rho Metropolis Hastings random walk std. deviation for rho
#' @param sigma2_diffuse magnitude defining "infinity" for diffuse initial condition
#' @param thin_rate thinning rate for posterior draws (every thin_rate sample is used)
#' @param log_state_only only log final state estimate (saves memory)
#' @export
sc_gibbs <- function(y,order_trend,order_cycle,freqs,mc_iter,
                     initial_sigma2_epsilon=1.e-4,
                     sigma2_zeta_fixed=FALSE,initial_sigma2_zeta=1.e-6,
                     initial_sigma2_kappa=1.e-11,sigma2_diffuse=1.e4,
                     rho_fixed=FALSE,initial_rho=0.9925,
                     rho_high=0.9945,
                     rho_low=0.99,
                     mh_sd_rho=0.000025,
                     thin_rate=5,
                     log_state_only=FALSE
                     )
{
  # Preliminaries constructing model and doing regression if it is offline
  nobs <- length(y) # total length of obs data, as T in the equation
  miss <- is.na(y)  # need this for computations on epsilon

  ygood <- y[!miss]
  ngood <- length(ygood)

  n_freq <- length(freqs)
  rho <- initial_rho
  scmat <- sc_model_mat(order_trend, order_cycle, freqs, rho)

  # V is observation error, 1x1
  V <- diag(c(1))

  # W is state innovation variance,
  Wtrend <- c(rep(0.,order_trend-1),1.)*initial_sigma2_zeta
  if (length(initial_sigma2_kappa) == 1){
    initial_sigma2_kappa <- rep(initial_sigma2_kappa,n_freq)
    Wcycle <- rep( c(rep(0,2*(order_cycle-1)),1,1),n_freq)*initial_sigma2_kappa
  }
  sigma2_kappa_itr <- initial_sigma2_kappa  # temporary holder during calcs
  W <- diag(c(Wtrend,Wcycle))

  nstate <- order_trend + 2*order_cycle*n_freq

  # calculate covariance matrix for nth-order cycle
  #covar_matrix_size = order_cycle**2
  #These are initial calculations that are updated as rho evolves
  covar_matrix_cycle <- list()
  covar_matrix_cycle_inv <- list()

  # These are to cache some covariance calculations for proposed new rho that might
  # not be accepted
  covar_new <- list()
  covar_inv_new <- list()

  for (jj in (1:n_freq)) {
    covar_matrix_cycle[[jj]] <- calc_cycle_covariance_matrix(order_cycle,rho,freqs[jj])
    covar_matrix_cycle_inv[[jj]] <- solve(covar_matrix_cycle[[jj]])
  }

  # Initial condition for the state is zero but evolves
  # todo: make sure this is true for all versions
  m0 <- rep(0, dim(W)[1])
  m0[1] <- 0.

  C0_trend <- diag(1,order_trend)*sigma2_diffuse
  # This tiles the covarinances for cycles diagonally
  # todo: scrutinize the initial_sigma2_kappa part
  C0_cycle <- dlm::bdiag(covar_matrix_cycle)*initial_sigma2_kappa
  C0 <- dlm::bdiag(C0_trend, C0_cycle)
  model <- dlm::dlm(FF=scmat$FF,V=V,GG=scmat$GG,W=W,m0=m0,C0=C0)

  # begin MCMC/Gibbs' sampling iteration
  print (paste('Start MCMC',date()))
  print (paste('Requested mc_iter =',mc_iter))
  flush.console()

  # These are for logging samples obtained during MCMC process
  state_thinned <- c()
  rho_save <- matrix(numeric(mc_iter),nrow=mc_iter,ncol=1)

  if (log_state_only){
    log_variances = FALSE
  }else{
    log_variances = TRUE
    #beta_draw_save <- matrix(numeric(mc_iter*k), nrow = mc_iter, ncol = k) # todo: removed intercept so k not k+1
    sigma2_epsilon_save <- numeric(mc_iter)
    sigma2_kappa_save <- matrix(numeric(mc_iter*n_freq), nrow = mc_iter, ncol = n_freq)
    sigma2_zeta_save <- numeric(mc_iter)
    # From here down also for logging, but assumes possible thinning
    subtide_thinned <- c()
    amplitudes_thinned <- list()
    phase_thinned <- list()
    freq_thinned <- list()
  }
  # Stores a tally of accepted/rejected M-H jumps
  accept_reject <- c(0,0)




  for (itr in(1:mc_iter)) {
    print("========")
    print(paste('MCMC iteration',itr))

    #step 1: observation equation.
    if (itr==1) {
      y_mod <- y - mean(y,na.rm=TRUE)  #todo: scrutinize need for this
    } else {
      state_part <-alpha_draw %*% t(model$FF)
      y_mod <- y - state_part
    }

    tide_res = y_mod[!miss]
    s2 <- crossprod(tide_res,tide_res) #equation (14.7) in BDA book
    dim(s2) <- NULL

    if(itr == 1){
      sigma2_eps_draw <- initial_sigma2_epsilon  #todo: hardwire
    }else{
       sigma2_eps_draw <- invgamma::rinvgamma(1,rate=s2/2.,shape=ngood/2.)
    }
    print(paste("New sigma2_epsilon ",sigma2_eps_draw))
    print(paste("sum squared residual",s2," s2/ngood=",s2/ngood))

    dlm::V(model) <- matrix(sigma2_eps_draw)
    if (log_variances){
      sigma2_epsilon_save[itr] <- sigma2_eps_draw
    }

    # Have to embed these back in the full size including missing data
    # The use of "res" is from an older variant with a linear predictor.
    # todo: remove if that capability is not added back
    res_good <- ygood
    res <- y # for size including missing data, don't care about y
    res [miss] <- NA
    res[!miss] <- res_good

    #step 2 - sample state (alpha) using dlmBSample
    modFilt <- dlm::dlmFilter(res, model)

    alpha_draw <- dlm::dlmBSample(modFilt)
    alpha_draw <- alpha_draw[-1,]
    m0_new <- alpha_draw[1,]
    state_part2 <-alpha_draw %*% t(model$FF)
    do_thinning <- (itr %% thin_rate) == 0
    if (do_thinning){    # todo: is this needed without the regression component?
      reconstruct <- state_part2
    }

    #step 3 - sample sigma_zeta (trend) and sigma_kappa (cycle),
    #         both follow inverse gamma distribution
    # notation in Harvey and Trimbur is beta(n) - beta(n-1) but we already used beta
    trend_diff <- alpha_draw[(2:nobs),order_trend]-alpha_draw[1:(nobs-1),order_trend]
    scale_zeta <- crossprod(trend_diff,trend_diff)
    shape_zeta <- nobs/2+.5

    if (sigma2_zeta_fixed){
      sigma2_zeta <- initial_sigma2_zeta
      print(paste("sigma2_zeta is fixed:",sigma2_zeta))
    }else{
      sigma2_zeta <- invgamma::rinvgamma(1, rate=scale_zeta/2., shape=shape_zeta)
      print(paste("New sigma2_zeta=",sigma2_zeta))
    }



    if (log_variances){
      sigma2_zeta_save[itr] <- sigma2_zeta
    }
    Wtrend <- c(rep(0.,order_trend-1),sigma2_zeta)
    print("*-*")

    #3.2 - cycles and kappa part
    Wcycle <- c()
    C0_cycle_new <- list()
    rho_old <- rho

    # Drawing rho regardless of whether it is fixed
    # wastes a few computer cycles, but drawing a rho
    # keeps the random seed consistent regardless of whether it is fixed.
    rho_new <- rho_old + rnorm(1,mean=0,sd=mh_sd_rho)

    rho_in_bounds <- ((rho_new >= rho_low) & (rho_new <= rho_high))
    if (rho_fixed){
      rho_in_bounds <- FALSE
      rho_new <- rho
    }

    pold <- 0
    pnew <- 0 #Accumulates log probability for rho
    sigma2_kappa_itr <- rep(0.,n_freq)
    for (jj in (1:n_freq)) {
      # each freq has its own variance
      freqs_sel <- freqs[jj]
      covar_matrix_inv_sel <- covar_matrix_cycle_inv[[jj]]
      # todo: this was rho_sel. double check this should be rho_old
      scale_kappa <- cal_shape_kappa(order_trend, order_cycle, jj, freqs_sel,
                                     rho_old, covar_matrix_inv_sel,
                                     nobs, alpha_draw)
      # todo: review scale
      # Draw sigma2_kappa
      shape_kappa <- nobs+order_cycle-1
      sigma2_kappa <- invgamma::rinvgamma(1, rate=scale_kappa/2., shape=shape_kappa)
      print(paste("freq jj =",jj," shape=",shape_kappa," scale=",scale_kappa,
                  " sigma2=",sigma2_kappa))

      if (log_variances){
        sigma2_kappa_save[itr,jj] <- sigma2_kappa
      }
      sigma2_kappa_itr[jj] <- sigma2_kappa
      Wcycle <- c(Wcycle,rep(0,2*(order_cycle-1)),sigma2_kappa,sigma2_kappa)

      if ((!rho_fixed) & rho_in_bounds){
        #
        covar_old<- covar_matrix_cycle[[jj]]
        covar_inv_old <- covar_matrix_cycle_inv[[jj]]
        arg_old <- scale_kappa

        pdivold <- 0.5*(determinant(sigma2_kappa*covar_old,
                                    logarithm=TRUE)$modulus[[1]])


        pold <- pold -0.5*arg_old-pdivold   # this is the log  probability


        covar_new[[jj]] <- calc_cycle_covariance_matrix(order_cycle,rho_new,freqs[jj])
        covar_inv_new[[jj]] <- solve(covar_new[[jj]])
        arg_new <- cal_shape_kappa(order_trend, order_cycle, jj, freqs_sel,
                                   rho_new, covar_inv_new[[jj]], nobs, alpha_draw)

        pdivnew <- 0.5*(determinant(sigma2_kappa*covar_new[[jj]],logarithm=TRUE)$modulus[[1]])
        pnew <- pnew -0.5*arg_new-pdivnew   # this is the log probability

      }
    }  # loop jj over frequencies

    #prior_old = dgamma(1.-rho_old,scale=0.000005,shape=600,log=TRUE)
    #prior_new = dgamma(1.-rho_new,scale=0.000005,shape=600,log=TRUE)
    pold <- pold  #+ prior_old
    pnew <- pnew  #+ prior_new

    # Now do the M-H step for rho
    rho <- rho_old    # This doesn't do anything, but it clarifies
    # If the new value is precluded by being out of bounds no action needed except to record it
    if ((!rho_fixed) & rho_in_bounds){
      if (pnew >= pold){
        # ratio > 1, accept without transformations and random draw
        print(paste("rho: accept ratio >= 1, accept new rho = ",rho_new))
        rho <- rho_new

      }else{
        logpdiff = pnew-pold
        accept_rate <- exp(logpdiff)
        u =runif(1,0.,1.)
        print(paste("rho: calc accept_rate=",accept_rate))
        if (u < accept_rate){
          print(paste("accept new rho = ",rho_new))
          rho <- rho_new
        }else{
          print(paste("rho reject based on ratio",rho_new," retain ",rho_old))
        }
      } # pnew >= pold
    }else{
      # If we get here, either we have fixed rho, it is rejected because it
      # is out of bounds (its prior is uniform) or it is rejected based
      # on posterior likelihood
      if (rho_fixed){
           rho <- rho_old
           print(paste("rho fixed, retain",rho_old))
        }else{
           print(paste("rho: rejected new rho:",rho_new," b/c out of bounds."," retain ",rho_old))
        }
    }

    if (rho != rho_old){
      if (rho !=  rho_new) { print("BUG in rho code")}
      rho_save[itr] <- rho_new
      accept_reject[1] <- accept_reject[1] + 1
      for (jj in 1:n_freq){
        covar_matrix_cycle[[jj]] <- covar_new[[jj]]
        covar_matrix_cycle_inv[[jj]] <- covar_inv_new[[jj]]
        C0_cycle_new[[jj]] <- covar_matrix_cycle[[jj]]*sigma2_kappa_itr[jj]
      }
    }else{
      rho <- rho_old # for clarity
      rho_save[itr] <- rho_old
      accept_reject[2] <- accept_reject[2] + 1
      for (jj in 1:n_freq){
        # todo: This is probably avoidable, but need to make sure C0_cycle_new is set for
        # the first iterations before a rho is accepted
        C0_cycle_new[[jj]] <- covar_matrix_cycle[[jj]]*sigma2_kappa_itr[jj]
      }
    }
    if (!rho_fixed){
      print("Accept/reject matrix:")
      print(accept_reject)
      print("*")
    }


    dlm::W(model) <- diag(c(Wtrend,Wcycle))
    C0_cycle <- dlm::bdiag(C0_cycle_new)
    dlm::C0(model) <- dlm::bdiag(C0_trend, C0_cycle)

    if (do_thinning){
      print(paste("iter=",itr))
      print ("thinned")
      state_thinned <- rbind(state_thinned,t(reconstruct))
      if (!log_state_only){
        subtide_thinned <-  rbind(subtide_thinned,t(alpha_draw[,1]))
        for (jj in (1:n_freq)) {
          flabels <- names(sc_freq)
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


  }

  print (paste('Finish MCMC',date()))
  flush.console()
  if (log_state_only){
    return(list(state_thinned=state_thinned))
  }else{
    return(list(sigma2_epsilon=sigma2_epsilon_save,
              sigma2_zeta=sigma2_zeta_save,
              sigma2_kappa=sigma2_kappa_save,
              rho = rho_save,
              order_trend=order_trend,
              order_cycle = order_cycle,
              freqs = freqs,
              state_thinned=state_thinned,
              subtide_thinned = subtide_thinned,
              amplitude_thinned=amplitudes_thinned,
              phase_thinned=phase_thinned,
              freq_thinned=freq_thinned))

    }
}
