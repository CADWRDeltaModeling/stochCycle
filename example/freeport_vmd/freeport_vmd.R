# This example has been tested against v09_m2n4_lmr_roughh_zetaf_kappndv
# It gives the same results for sigma2_zeta_fixed = FALSE when the prior is zero/improper
# but the option isn't available in v09.

library(stochCycle)
require(VMDecomp)
library(pracma)
library(seewave)


CFS2CMS = 0.028316847
# The data are loaded in a unit that is large (cubic feet per second)
# This transformation to makes it nearly unit scaled
fpt = freeport_discharge*CFS2CMS/100.

# list of tidal frequencies in radians per hour
samples_hr = 4. # data are 15min
tide_freqs = tide_freq_radhr/samples_hr

# this subsets approximately six months, minimal missing data
# create some missing data
select <- 5000:23000
y <- fpt[select,1]



default_vmd_params <- list(alpha = 2000,
                           tau = 0,
                           DC = FALSE,
                           init = 1,
                           tol = 1e-6)

#res_k <- estimate_k_modes(signal_1d = y,
#                          cor_thresh = 0.01,
#                          default_vmd_params = default_vmd_params,
#                          min_K = 2,
#                          seed = 1,
#                          verbose = TRUE)

K <- 8 # res_k
res_1d <- vmd(data = y,
              alpha = 4000,
              tau = 0,
              K = K,
              DC = FALSE,
              init = 1,
              tol = 1e-6,
              verbose = FALSE)

#plotting intrinsic mode functions (IMF)
op <- par(mfrow = c(ceiling(K/2), 2))
for (item in 1:5) {
  item_mode = glue::glue("IMF {item}")
  plot(x = res_1d$u[1:10000, item], type = 'l', main = item_mode, xlab = 'Time', ylab = '')
}

#spectral analysis
op <- par(mfrow = c(ceiling(K/3), 3))
sampling_interval_day <- 15.0/(24*60) # sampling duration is 15 minutes
for (item in 1:K) {
  x.spec <- spectrum(res_1d$u[,item], log="no", span=10,plot=FALSE)
  spx <- x.spec$freq/sampling_interval_day
  spy <- 2*x.spec$spec
  item_mode = glue::glue("Mode {item}")
  plot(spy~spx,xlab="frequency (cpd)",main = item_mode, ylab="spectral density",type="l", xlim=c(0,10))
}

#central frequency
central_freq <- tail(res_1d$omega,1)/sampling_interval_day
central_period <- 1/central_freq
print ("central frequency (cpd)")
print (central_freq)
print ("central period (day)")
print (central_period)

#calculate mean amplitudes and periods of VMD modes
#according to Gan et al. 2021 paper
op <- par(mfrow = c(ceiling(K/2), 2))
total_period_day <- (length(y)-1)*15/(1440.0)
sample_frq <- 1/(15*60) #as Hz
summary <- data.frame()
envelop <- 1:length(y)/4

last_K <- 5 # based on analysis, these are the interpretable ones
for (item in 1:last_K){
  mode <- item
  item_mode = glue::glue("Mode {item}")
  temp <- findpeaks(res_1d$u[,item])
  no_peak <- dim(temp)[1]
  mean_period_day <- total_period_day/no_peak
  mean_amp <- mean(temp[,1])
  new_row <- data.frame(mode, no_peak,mean_period_day, mean_amp)
  summary <- rbind(summary, new_row)
  hilb <- hilbert(res_1d$u[,item],sample_frq)
  amp <- as.vector(abs(hilb))
  phase <- as.vector(angle(hilb))
  if (item == 1){
    newcol <- amp
  }else{
    newcol <- cbind(amp,phase)
  }
  newcols <- data.frame(newcol)

  envelop <- cbind(envelop, newcols)
  plot(phase,type='l',main = item_mode, xlab = 'Time', ylab = 'Amplitude')
  points(temp[,2],temp[,1])
}
envelop[,2] <- res_1d$u[,1] # subtide not an envelop
names(envelop) <- c('time','subtide',
                    'D2_amp','D2_phase',
                    'D3_amp','D3_phase',
                    'D4_amp','D4_phase',
                    'D1_amp','D1_phase')
print("")
print("Mean amplitudes and periods of VMD modes")
print(summary)
write.csv(format(envelop,digits=3,scientific=FALSE) ,file="freeport_vmd.csv",quote=FALSE,row.names=FALSE)

ts.plot(envelop[1:500,'D2_phase'])


par(mfrow=c(1,1))
ts.plot(y)
lines(envelop[,'subtide'],col='red')


