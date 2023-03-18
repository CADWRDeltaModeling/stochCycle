library(stochCycle)
require(VMDecomp)
library(pracma)
library(seewave)

test_input <- read.csv("../test_tide.csv",header=FALSE,stringsAsFactors = FALSE)[,2]
test <- test_input

# Note that y already has random component from the python script with SD = 0.01
y <- test


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

K <- 5 # res_k
res_1d <- vmd(data = y,
              alpha = 2000,
              tau = 0,
              K = K,
              DC = FALSE,
              init = 1,
              tol = 1e-6,
              verbose = FALSE)

#plotting intrinsic mode functions (IMF)
op <- par(mfrow = c(ceiling(K/2), 2))
for (item in 1:K) {
  item_mode = glue::glue("IMF {item}")
  plot(x = res_1d$u[, item], type = 'l', main = item_mode, xlab = 'Time', ylab = '')
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
envelop <- 1:length(y)
for (item in 1:K) {
  mode <- item
  item_mode = glue::glue("Mode {item}")
  temp <- findpeaks(res_1d$u[,item])
  no_peak <- dim(temp)[1]
  mean_period_day <- total_period_day/no_peak
  mean_amp <- mean(temp[,1])
  new_row <- data.frame(mode, no_peak,mean_period_day, mean_amp)
  summary <- rbind(summary, new_row)
  hilb <- hilbert(res_1d$u[,item],sample_frq)
  env <- as.vector(abs(hilb))
  new_col <- data.frame(env)
  envelop <- cbind(envelop, new_col)
  plot(env,type='l',main = item_mode, xlab = 'Time', ylab = 'Amplitude')
  points(temp[,2],temp[,1])
}
res<-1:length(y)
sub <- envelop[,6]
res <- cbind(res, sub)
D1_amp <- envelop[,4]
res <- cbind(res, D1_amp)
D2_amp <- envelop[,1]
res <- cbind(res, D2_amp)
D3_amp <- envelop[,5]
res <- cbind(res, D3_amp)
D4_amp <- envelop[,3]
res <- cbind(res, D4_amp)

names(res) <- c('time','subtide','D1','D2','D3','D4')
print("")
print("Mean amplitudes and periods of VMD modes")
print(summary)
write.csv(format(res,digits=3,scientific=FALSE) ,file="tropical_vmd_amp.csv",quote=FALSE,row.names=FALSE)




