# This example has been tested against v09_m2n4_lmr_roughh_zetaf_kappndv
# It gives the same results for sigma2_zeta_fixed = FALSE when the prior is zero/improper
# but the option isn't available in v09.

library(stochCycle)
#library(dlm)
#library(invgamma)
#library(MASS)
#library(lqmm)

CFS2CMS <- 0.028316847
# The data are loaded in a unit that is large (cubic feet per second)
# This transformation to makes it nearly unit scaled
fpt_input <- read.csv("../test_tide.csv",header=FALSE,stringsAsFactors = FALSE)[,2]
fpt <- fpt_input  #*CFS2CMS/100.

# list of tidal frequencies in radians per hour
samples_hr = 4. # data are 15min
tide_freqs = tide_freq_radhr/samples_hr

y <- fpt #[,1]


regress_names <- c("K1","M2","O1","S2","Q1","N2","L2","M4","MK3","MO3")
regress_freq <- tide_freqs[regress_names]

# Others are M05(one cycle per two tidal days, inertial
sc_names <- c("M1","M2","M3","M4")
sc_freq <- tide_freqs[sc_names]

mc_iter <- 2001
regress_names <- c("K1","M2","O1","S2","Q1","N2","L2","M4","MK3","MO3") #,"M6","MK5")
order_trend <- 2
order_cycle <- 3

# Leftover
zeta_prior_shape=0.
zeta_prior_rate <- 0.


# This enables an alternate prior that is concentrated at high values
# It isn't parameterized yet
use_rho_prior <- 0.
regress <- FALSE
regress_lmr <- FALSE
sigma2_zeta_fixed <- FALSE   # previously not fixed
#1 0.001 works for flow(3,3), becomes smaller 0.1 and (2,2)
initial_sigma2_zeta <- 2.5e-07   # from successful 2-4 field example
initial_sigma2_kappa <- 1.e-13 # previously 1.35e-11 for 2-4
initial_sigma2_epsilon <- 1.e-3  #was e-6


# For 2-3
mh_sd_rho <- 0.0002
mh_sd_kappa <- 5.e-13
mh_sd_zeta <- 1.e-8
mh_sd_eps <- 0.0001 # was 1e-6 not -5
initial_rho <- 0.99
initial_sigma2_zeta <- 2.e-07   # from successful 2-4 field example
initial_sigma2_kappa <- 1.35e-12 # previously 1.35e-11 for 2-4

all_test_freq <- regular_filter_targets_23()
test_freqs <-  all_test_freq$test_freqs
test_targs <- all_test_freq$test_targets
filter_scale <- 0.01 # 0.05 was good for 2-4 but not for 2-3


sigma2_diffuse <- 100. #
thin_rate <- 5
max_rho <- 0.9975

sc_out <- sc_gibbs_filter_prior(y,order_trend,order_cycle,
          sc_freq,test_freqs,test_targs,filter_scale,samples_hr,
          regress,regress_lmr,regress_freq,initial_sigma2_epsilon,
          use_rho_prior,initial_rho,mh_sd_rho,max_rho,
          sigma2_zeta_fixed,zeta_prior_rate,zeta_prior_shape,mh_sd_zeta,mh_sd_eps,
          initial_sigma2_zeta,initial_sigma2_kappa,mh_sd_kappa,sigma2_diffuse,
          mc_iter,thin_rate)

mod <- sc_model_from_gibbs(sc_out)
modFilt <- dlm::dlmFilter(y, mod)
modSmooth <- dlm::dlmSmooth(modFilt)
recon <- reconstruct_state(modSmooth$s,order_trend,order_cycle,4)
do_plot(ysave,y,recon,500:2100)
interp <- interpret_state(modSmooth$s,order_trend,order_cycle,length(sc_freq),
                          flabels=sc_names)


param_interp_amp<-cbind(interp$subtide,interp$amplitude$M1,interp$amplitude$M3)
ts.plot(param_interp_amp,ylim=c(-0.25,0.9))
write.csv(x=param_interp_amp,file="SCHA_2_3_param_amp.csv")

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

