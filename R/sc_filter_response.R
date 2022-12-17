

lowpass_gain_term <- function(lambdaf,m,q){
  elambdaf <- -1.i*lambdaf
  elambdaf <- exp(elambdaf)
  gain0 <- `^`((1.-elambdaf)*(1.-1./elambdaf),m)
  gain1 <- Re(q/gain0)
}

spectral_density_real <- function(lambdaf,lambdac,rho,nn,q){
  n = as.double(nn)
  rho2 = rho*rho
  num = lambdaf*0.
  for (j in 0:(nn)){
    for (k in 0:(nn)){
      jkdiff <- as.double(j-k)
      jksum <- as.double(j+k)
      term <- ((-1.)^jksum)*choose(nn,j)*choose(nn,k)*cos(lambdac*jkdiff)*cos(lambdaf*jkdiff)*(rho^jksum)
      num <- num + term
    }
  }

  denom= 1.+ 4.*rho2*cos(lambdac)**2. + rho2*rho2 - 4.*rho*(1.+rho2)*cos(lambdac)*cos(lambdaf)+2.*rho2*cos(2.*lambdaf)
  denom = denom**n
  two_pi = pi*2.
  if(q > 0){
    # return power spectrum (Trimbur 2006 pg 11, no eq number)
    g = (q/two_pi)*(num/denom)
    return(g)
  }


  g_without_sigma_term = num/(two_pi*denom)

  #return spectral density
  # these are the numerator and denominator of the variance (16) ignoring
  # q or sigma_kappa because it divides out
  # todo: need to check the range of the sum and the power rho_to_2i as there were inconsistencies
  var_den = (1.-rho2)**(2.*nn-1.)
  var_num = 0.
  for (i in 0:(nn-1)){    # todo is the -1 right?
    rho_to_2i = rho2**as.double(i)
    var_num <- var_num + (choose(nn-1,i)*rho_to_2i)**2.
  }
  var_without_sigma_term = var_num/var_den
  return(g_without_sigma_term/var_without_sigma_term)
}

#' Evaluate filter gain at requested frequencies
#'
#' @param evalfreqs vector of frequencies where gain is requested
#' @param freqc list of frequencies for cycles
#' @param m order of trend
#' @param n order of cycle
#' @param rho dampening
#' @param qzeta variance of trend divided by variance of observation
#' @param qkappa variance of cycle divided by variance of observation
#' @return vector of gains
#' @export
calculate_gain<-function(evalfreqs,freqc,m,n,rho,qzeta,qkappa){
  lowgain<-lowpass_gain_term(evalfreqs,m,qzeta)
  nfreq <- length(freqc)
  if (is.atomic(qkappa)){qkappa<-rep(qkappa,nfreq)}
  gain <- list()
  for (i in 1:nfreq){
    freq = freqc[[i]]
    gain[[i]] <- spectral_density_real(evalfreqs,freq,rho,n,qkappa[i])
  }
  gain[[nfreq+1]]<-lowgain
  return(gain)
}

#' Plot gains of the sc components
#'
#' @param sc_freq list of frequencies for cycles
#' @param m order of trend
#' @param n order of cycle
#' @param rho dampening
#' @param qzeta variance of trend divided by variance of observation
#' @param qkappa variance of cycle divided by variance of observation
#' @param nsample number of evaluaiton points up to 5cpd #todo: flexible
#' @return vector of gains
#' @export
plot_gain<- function(sc_freq,m,n,rho,qzeta,qkappa,nsample){

  cycles_per_day <- seq(1/nsample,5,5/nsample)
  rad_hr <- 2*pi*cycles_per_day/24.
  nfreq <- length(sc_freq)
  if (is.atomic(qkappa)){qkappa<-rep(qkappa,nfreq)}

  sd <- calculate_gain(rad_hr,sc_freq,m,n,rho,qzeta,qkappa)
  nall <- length(sd)

  total <- Reduce('+',sd)


  total <- total + 1.
  for (i in 1:nall){
    sd[[i]] <- sd[[i]]/total
  }

  plot(cycles_per_day,sd[[1]],ylim=c(0,1),xlim=c(0.,4),type='l',col='green')
  lines(cycles_per_day,sd[[2]],col='darkblue')
  lines(cycles_per_day,sd[[3]],col='maroon')
  lines(cycles_per_day,sd[[4]],col='purple')
  lines(cycles_per_day,sd[[5]],col='red')
  return(sd)
}


#' Evaluate improper prior representing filter quality
#'
#' @param test_freqs list of frequencies at which filter quality is to be eval'd
#' @param test_targets ideal gain at test freqs (1 or 0)
#' @param scale scale of penalty. The smaller, the less forgiving
#' @param sc_freq list of frequencies for cycles
#' @param m order of trend
#' @param n order of cycle
#' @param rho dampening
#' @param qzeta variance of trend divided by variance of observation
#' @param qkappa variance of cycle divided by variance of observation
#' @return cost function evaluated
#' @export
filter_prior <- function(test_freqs,test_targets,scale,
                         sc_freq,m,n,rho,qzeta,qkappa){
  nfreq <- length(sc_freq)
  if (is.atomic(qkappa)){qkappa<-rep(qkappa,nfreq)}
  #print(paste("cost eval scale=",scale,"rho=",rho,"qzeta=",qzeta,
  #            "qkappa=",qkappa[1],qkappa[2],qkappa[3],qkappa[4],
  #             "number of test freqs:",length(test_freqs)))
  cost <- 0.
  for (ii in 1:length(test_freqs)){
    tester <- test_freqs[[ii]]
    target <- test_targets[[ii]]
    trad_per_hr <- 2*pi*tester/24.
    sd <- calculate_gain(trad_per_hr,sc_freq,m,n,rho,qzeta,qkappa)
    nall <- length(sd)

    total <- Reduce('+',sd)
    total <- total + 1.
    sd[[ii]] <- sd[[ii]]/total
    ntarget <- length(target)
    for (jjj in 1:ntarget){
      if (target[jjj] == 1.){
        sd[[ii]][jjj] <- 1. - sd[[ii]][jjj]
      }
      add_cost <- (sd[[ii]][jjj]/scale)**2.
      cost <- cost + add_cost

    }
  }
  #print(paste("Return cost:",format(cost,scientific=TRUE)))

  return(cost)
}

#' Return a list of targets and frequencies that is good for D1-D2-D3-D4
#' @return list of testing frequencies and ideal filter targets (0 or 1)
#' @export
sparse_filter_targets <- function(){
  # test frequencies in cpd for evaluation of log prior based on filter performance
  testfreqs <- list()
  test_targets <- list()
  testfreqs[[1]] <- c(0.4,0.45,0.5,0.525,0.55,0.9,0.93,0.96,1.0,1.04,1.08,1.5)
  test_targets[[1]] <- c(0.,0.,0.,0.,0.,1.,1.,1.,1.,1.,1.0,0.)
  # 1.932 is Mw
  testfreqs[[2]] <- c(0.5,1.,1.5,1.8,2.,2.05,2.5)
  test_targets[[2]] <- c(0.,0.,0.,1.,1.,1.,0.)
  testfreqs[[3]] <- c(0.5,1.,1.5,1.95,2.5,2.8,3.,3.05,3.5)
  test_targets[[3]] <- c(0.,0.,0.,0.,0.,1.,1.,1.,0.)
  testfreqs[[4]] <- c(0.5,1.,1.95,3.,3.7,4.,4.1,4.5)
  test_targets[[4]] <- c(0.,0.,0.,0.,1.,1.,1.,0.)
  testfreqs[[5]] <- c(0.1,0.2,0.3,0.4,0.45,0.5,0.525,0.9,0.95,1.0,1.05,1.95,2.)
  test_targets[[5]] <- c(1.,1.,1.,1.,1.,1.,1.,0.,0.,0.,0.,0.)
  out <- list(test_freqs=testfreqs,test_targets=test_targets)
}

#' Return a list of targets and frequencies that is good for D1-D2-D3-D4
#' @return list of testing frequencies and ideal filter targets (0 or 1)
#' @export
regular_filter_targets <- function(){
  # test frequencies in cpd for evaluation of log prior based on filter performance
  testfreqs <- list()
  test_targets <- list()
  testfreqs[[1]] <- c(0.4,0.45,0.5,0.525,0.55,0.9,0.93,0.96,1.0,1.04,1.08,1.5)
  test_targets[[1]] <- c(0.,0.,0.,0.,0.,1.,1.,1.,1.,1.,1.,0.)
  testfreqs[[2]] <- c(0.5,1.,1.5,1.85,2.,2.15,2.5)
  test_targets[[2]] <- c(0.,0.,0.,1.,1.,1.,0.)
  testfreqs[[3]] <- c(0.5,1.,1.5,2.,2.5,2.85,3.,3.15,3.5)
  test_targets[[3]] <- c(0.,0.,0.,0.,0.,1.,1.,1.,0.)
  testfreqs[[4]] <- c(0.5,1.,2.,3.,3.8,4.,4.15)
  test_targets[[4]] <- c(0.,0.,0.,0.,1.,1.,1.)
  testfreqs[[5]] <- c(0.1,0.2,0.3,0.4,0.45,0.5,0.525,0.55,0.9,0.93,0.96,1.0,1.04,1.08,1.5,2.,2.5)
  test_targets[[5]] <- c(1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.)
  out <- list(test_freqs=testfreqs,test_targets=test_targets)
}

#' Return a list of targets and frequencies that is good for D1-D2-D3-D4
#' @return list of testing frequencies and ideal filter targets (0 or 1)
#' @export
regular_filter_targets_more1 <- function(){
  # test frequencies in cpd for evaluation of log prior based on filter performance
  testfreqs <- list()
  test_targets <- list()
  testfreqs[[1]] <- c(0.35,0.4,0.45,0.5,0.9,0.92,0.94,0.96,0,98,1.0,1.02,1.04,1.06,1.08,1.5)
  test_targets[[1]] <- c(0.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.)
  testfreqs[[2]] <- c(0.5,1.,1.5,1.85,2.,2.15,2.5)
  test_targets[[2]] <- c(0.,0.,0.,1.,1.,1.,0.)
  testfreqs[[3]] <- c(0.5,1.,1.5,2.,2.5,2.85,3.,3.15,3.5)
  test_targets[[3]] <- c(0.,0.,0.,0.,0.,1.,1.,1.,0.)
  testfreqs[[4]] <- c(0.5,1.,2.,3.,3.8,4.,4.15)
  test_targets[[4]] <- c(0.,0.,0.,0.,1.,1.,1.)
  testfreqs[[5]] <- c(0.1,0.2,0.3,0.35,0.9,0.93,0.96,1.0,1.04,1.08,1.85,1.95)
  test_targets[[5]] <- c(1.,1.,1.,1.,0.,0.,0.,0.,0.,0.,0.,0.)
  out <- list(test_freqs=testfreqs,test_targets=test_targets)
}

#' Return a list of targets and frequencies that is good for D1-D2-D3-D4
#' @return list of testing frequencies and ideal filter targets (0 or 1)
#' @export
regular_filter_targets_more1b <- function(){
  # test frequencies in cpd for evaluation of log prior based on filter performance
  testfreqs <- list()
  test_targets <- list()
  testfreqs[[1]] <- c(0.35,0.4,0.45,0.5,0.9,0.92,0.94,0.96,0,98,1.0,1.02,1.04,1.06,1.08,1.5)
  test_targets[[1]] <- c(0.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.)
  testfreqs[[2]] <- c(0.5,1.,1.5,1.85,2.,2.15,2.5)
  test_targets[[2]] <- c(0.,0.,0.,1.,1.,1.,0.)
  testfreqs[[3]] <- c(0.5,1.,1.5,2.,2.5,2.85,3.,3.15,3.5)
  test_targets[[3]] <- c(0.,0.,0.,0.,0.,1.,1.,1.,0.)
  testfreqs[[4]] <- c(0.5,1.,2.,3.,3.8,4.,4.15)
  test_targets[[4]] <- c(0.,0.,0.,0.,1.,1.,1.)
  testfreqs[[5]] <- c(0.1,0.2,0.3,0.35,0.93,0.96,1.0,1.04,1.08,1.45,1.5,1.55,1.85,1.95)
  test_targets[[5]] <- c(1.,1.,1.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.)
  out <- list(test_freqs=testfreqs,test_targets=test_targets)
}


#' Return a list of targets and frequencies that is good for D1-D2-D3-D4
#' @return list of testing frequencies and ideal filter targets (0 or 1)
#' @export
regular_filter_targets_23 <- function(){
  # test frequencies in cpd for evaluation of log prior based on filter performance
  testfreqs <- list()
  test_targets <- list()
  testfreqs[[1]] <- c(0.35,0.4,0.45,0.96,0,98,1.0,1.02,1.04,1.06,1.08,1.5)
  test_targets[[1]] <- c(0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,0.)
  testfreqs[[2]] <- c(0.5,1.,1.5,1.85,2.,2.15,2.5)
  test_targets[[2]] <- c(0.,0.,0.,1.,1.,1.,0.)
  testfreqs[[3]] <- c(0.5,1.,1.5,2.,2.5,2.85,3.,3.15,3.5)
  test_targets[[3]] <- c(0.,0.,0.,0.,0.,1.,1.,1.,0.)
  testfreqs[[4]] <- c(0.5,1.,2.,3.,3.8,4.,4.15)
  test_targets[[4]] <- c(0.,0.,0.,0.,1.,1.,1.)
  testfreqs[[5]] <- c(0.1,0.2,0.3,0.96,0.98,1.0,1.04,1.08,1.45,1.5,1.55,1.85,1.95,2.5)
  test_targets[[5]] <- c(1.,1.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.)
  out <- list(test_freqs=testfreqs,test_targets=test_targets)
}

#' Return a list of targets and frequencies that is good for D1-D2-D3-D4
#' @return list of testing frequencies and ideal filter targets (0 or 1)
#' @export
regular_filter_targets_23b <- function(){
  # test frequencies in cpd for evaluation of log prior based on filter performance
  testfreqs <- list()
  test_targets <- list()
  testfreqs[[1]] <- c(0.35,0.4,0.45,0.96,0,98,1.0,1.02,1.04,1.06,1.08,1.5)
  test_targets[[1]] <- c(0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,0.)
  testfreqs[[2]] <- c(0.5,1.,1.5,1.85,2.,2.15,2.5)
  test_targets[[2]] <- c(0.,0.,0.,1.,1.,1.,0.)
  testfreqs[[3]] <- c(0.5,1.,1.5,2.,2.5,2.85,3.,3.15,3.5)
  test_targets[[3]] <- c(0.,0.,0.,0.,0.,1.,1.,1.,0.)
  testfreqs[[4]] <- c(0.5,1.,2.,3.,3.8,4.,4.15)
  test_targets[[4]] <- c(0.,0.,0.,0.,1.,1.,1.)
  testfreqs[[5]] <- c(0.1,0.2,0.3,0.96,0.98,1.0,1.04,1.08,1.45,1.5,1.55,1.85,1.95,2.5,7.5)
  test_targets[[5]] <- c(1.,1.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.)
  out <- list(test_freqs=testfreqs,test_targets=test_targets)
}

#' Return a list of targets and frequencies that is good for D1-D2-D3-D4
#' @return list of testing frequencies and ideal filter targets (0 or 1)
#' @export
regular_filter_targets_24b <- function(){
  # test frequencies in cpd for evaluation of log prior based on filter performance
  testfreqs <- list()
  test_targets <- list()
  testfreqs[[1]] <- c(0.35,0.4,0.45,0.96,0,98,1.0,1.02,1.04,1.06,1.08,1.5)
  test_targets[[1]] <- c(0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,0.)
  testfreqs[[2]] <- c(0.5,1.,1.5,1.85,2.,2.15,2.5)
  test_targets[[2]] <- c(0.,0.,0.,1.,1.,1.,0.)
  testfreqs[[3]] <- c(0.5,1.,1.5,2.,2.5,2.85,3.,3.15,3.5)
  test_targets[[3]] <- c(0.,0.,0.,0.,0.,1.,1.,1.,0.)
  testfreqs[[4]] <- c(0.5,1.,2.,3.,3.8,4.,4.15)
  test_targets[[4]] <- c(0.,0.,0.,0.,1.,1.,1.)
  testfreqs[[5]] <- c(0.1,0.2,0.3,0.4,0.94,0.96,0.98,1.0,1.04,1.08,1.45,1.5,1.55,1.85,1.95,2.5,3.5,4.,4.5,5.0,5.5,6.0,6.5,7.,7.5,8.)
  test_targets[[5]] <- c(1.,1.,1.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.0,0.0,0.0,0.,0.,0.,0.)
  out <- list(test_freqs=testfreqs,test_targets=test_targets)
}

#' Return a list of targets and frequencies that is good for D1-D2-D3-D4
#' @return list of testing frequencies and ideal filter targets (0 or 1)
#' @export
regular_filter_targets_1thru6 <- function(){
  # test frequencies in cpd for evaluation of log prior based on filter performance
  testfreqs <- list()
  test_targets <- list()
  testfreqs[[1]] <-    c(0.35,0.4,0.45,0.96,0,98,1.0,1.02,1.04,1.06,1.08,1.5,2.0)
  test_targets[[1]] <- c(0.00,0.0,0.00,1.00,1.00,1.0,1.00,1.00,1.00,1.00,0.0,0.0)
  testfreqs[[2]] <-    c(0.5,1.,1.5,1.85,2.,2.15,2.5,3.,3.5,4.0,4.5)
  test_targets[[2]] <- c(0.,0.,0.,1.,1.,1.,0.,0.,0.,0.,0.)
  testfreqs[[3]] <-    c(0.5,1.,1.5,2.,2.5,2.85,3.,3.15,3.5,4.0)
  test_targets[[3]] <- c(0.0,0.,0.0,0.,0.0,1.00,1.,1.00,0,0,0.0)
  testfreqs[[4]]    <- c(0.5,1.,2.0,3.0,3.8,4.0,4.15,4.5)
  test_targets[[4]] <- c(0.0,0.,0.0,0.0,1.0,1.0,1.00,0.0)
  testfreqs[[5]] <-    c(0.5,1.,2.,3.,4.,4.5,5.)
  test_targets[[5]] <- c(0.0,0.,0.,0.,0.,0.0,1.)
  testfreqs[[6]] <-    c(0.5,1.,2.,3.,4.,5.,6.)
  test_targets[[6]] <- c(0.0,0.,0.,0.,0.,0.,1.)
  testfreqs[[7]] <-    c(0.1,0.2,0.3,0.4,0.94,0.96,0.98,1.0,1.04,1.08,1.45,1.5,1.55,1.85,1.95,2.5,3.5,4.,4.5,5.0,5.5,6.0,6.5,7.,7.5,8.)
  test_targets[[7]] <- c(1.0,1.0,1.0,1.0,0.00,0.00,0.00,0.0,0.00,0.00,0.00,0.0,0.00,0.00,0.00,0.0,0.0,0.,0.0,0.0,0.0,0.0,0.0,0.,0.0,0.)
  out <- list(test_freqs=testfreqs,test_targets=test_targets)
}


#' Example of filter prior/loss function evaluation
#' @return cost function evaluated
#' @export
example_filter <- function(eps,zeta,kappa,rho,m=2,n=4,scale=0.001,nsample=2400){
  tide_freqs = tide_freq_radhr
  sc_names <- c("M1","M2","M3","M4")
  sc_freq <- tide_freqs[sc_names]
  sigma2_zeta <- zeta
  sigma2_eps <- eps
  sigma2_kappa <- kappa
  qzeta <- sigma2_zeta/sigma2_eps
  qkappa <- sigma2_kappa/sigma2_eps

  # test frequencies for evaluation of log prior based on filter performance
  alltargs <-regular_filter_targets()
  ttargets <- alltargs$test_targets
  testfreqs <- alltargs$test_freqs

  plot_gain(sc_freq,m,n,rho,qzeta,qkappa,nsample)
  cost <-filter_prior(testfreqs,ttargets,scale,
                      sc_freq,m,n,rho,qzeta,qkappa)


  print(paste("Cost:",cost))
  format(cost,scientific=TRUE)
}

