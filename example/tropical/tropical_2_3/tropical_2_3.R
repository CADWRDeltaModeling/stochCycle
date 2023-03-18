
library(stochCycle)

CFS2CMS <- 0.028316847
# The data are loaded in a unit that is large (cubic feet per second)
# This transformation to makes it nearly unit scaled
test_input <- read.csv("../test_tide.csv",header=FALSE,stringsAsFactors = FALSE)[,2]
test <- test_input  #*CFS2CMS/100.

# list of tidal frequencies in radians per hour
samples_hr = 4. # data are 15min
tide_freqs = tide_freq_radhr/samples_hr

# Note that y already has random component from the python script with SD = 0.01
y <- test


# Others are M05(one cycle per two tidal days, inertial
sc_names <- c("M1","M2","M3","M4")
sc_freq <- tide_freqs[sc_names]

mc_iter <- 2001
order_trend <- 2
order_cycle <- 3

rho_fixed <- TRUE
initial_rho <- 0.997
sigma2_zeta_fixed <- FALSE
initial_sigma2_zeta=1.e-10
initial_sigma2_kappa=1.35e-13

sc_out <- sc_gibbs(y,order_trend,order_cycle,sc_freq,mc_iter,
                   rho_fixed = rho_fixed,initial_rho=initial_rho,
                   initial_sigma2_epsilon=0.004,
                   sigma2_zeta_fixed=sigma2_zeta_fixed,
                   initial_sigma2_zeta=initial_sigma2_zeta,
                   initial_sigma2_kappa)


mod <- sc_model_from_gibbs(sc_out)
modFilt <- dlm::dlmFilter(y, mod)
modSmooth <- dlm::dlmSmooth(modFilt)
recon <- reconstruct_state(modSmooth$s,order_trend,order_cycle,4)
do_plot(y,y,recon,500:2100)
interp <- interpret_state(modSmooth$s,order_trend,order_cycle,length(sc_freq),
                          flabels=sc_names)

param_interp_amp<-cbind(interp$subtide,interp$amplitude$M1,interp$amplitude$M3)
ts.plot(param_interp_amp,ylim=c(-2,2))
write.csv(x=param_interp_amp,file="SCHA_2_4_param_amp.csv")

ts.plot(sc_out$amplitude_thinned$M2[6,],ylim=c(0.,1.5))
lines(sc_out$amplitude_thinned$M1[6,],col="blue")
lines(sc_out$amplitude_thinned$M3[6,],col="red")
lines(sc_out$amplitude_thinned$M4[22,],col="purple")

D1anal <- aggregate_iter(sc_out$amplitude_thinned$M1,nburn=-1)
ts.plot(D1anal,col="blue",ylim=c(0.,1.5))
D2anal <- aggregate_iter(sc_out$amplitude_thinned$M2,nburn=-1)
lines(D2anal,col="black")
D3anal <- aggregate_iter(sc_out$amplitude_thinned$M3,nburn=-1)
lines(D3anal,col="red")
D4anal <- aggregate_iter(sc_out$amplitude_thinned$M4,nburn=-1)
lines(D4anal,col="green")

ts.plot(recon$subtide[400:2100])
length(recon$subtide)
D1anal <- aggregate_iter(sc_out$amplitude_thinned$M1,nburn=-1)
ts.plot(D1anal[500:2100],col="blue",ylim=c(0.,1.))
lines(interp$amplitude$M1[500:2100])
