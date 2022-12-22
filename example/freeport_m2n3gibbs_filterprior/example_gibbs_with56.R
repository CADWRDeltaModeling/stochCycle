# This example has been tested against v09_m2n4_lmr_roughh_zetaf_kappndv
# It gives the same results for sigma2_zeta_fixed = FALSE when the prior is zero/improper
# but the option isn't available in v09.

library(stochCycle)

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
y <- fpt[,1]
ysave <- fpt[select,]
y <- y[select]



regress_names <- c("K1","M2","O1","S2","Q1","N2","L2","M4","MK3","MO3")
regress_freq <- tide_freqs[regress_names]

# Others are M05(one cycle per two tidal days, inertial
sc_names <- c("M1","M2","M3","M4","M5","M6")
sc_freq <- tide_freqs[sc_names]

mc_iter <- 1001
regress_names <- c("K1","M2","O1","S2","Q1","N2","L2","M4","MK3","MO3") #,"M6","MK5")
order_trend <- 2
order_cycle <- 3

# Leftover
zeta_prior_shape=0.
zeta_prior_rate <- 0.

use_rho_prior <- 0.
regress <- FALSE
regress_lmr <- FALSE
sigma2_zeta_fixed <- FALSE   # previously not fixed
sigma2_diffuse <- 100. #
thin_rate <- 5
max_rho <- 0.996



mh_sd_rho <- 0.0002
mh_sd_kappa <- 1.e-13
mh_sd_zeta <- 1.e-8
mh_sd_eps <- 0.00001
initial_rho <- 0.995
initial_sigma2_zeta <- 1.e-08   # from successful 2-4 field example
initial_sigma2_kappa <- 1.35e-13 # previously 1.35e-11 for 2-4
initial_sigma2_epsilon <- 5.e-4


# Restart
initial_sigma2_epsilon <- 0.000432
initial_sigma2_zeta <- 2.1e-07
inital_rho <- 0.9959
mh_sd_rho <- 0.00025
if (TRUE){
  initial_sigma2_kappa <- c(9.24e-12,1.697e-12,
                            4.982e-12,2.04e-09,
                            1.45e-13,1.47e-13)

  mh_sd_kappa <- c(5e-13,1e-13,1e-13,5e-11,5e-15,5e-15)
  mh_sd_zeta <- 1.e-8
  mh_sd_eps <- 0.00001

}

# This contains test frequencies for species 1 through 6
all_test_freq <- regular_filter_targets_1thru6()
test_freqs <-  all_test_freq$test_freqs
test_targs <- all_test_freq$test_targets
filter_scale <- 0.004  # was 0.02 in run with D8



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


# Collect output
ncycle <- length(sc_freq)
recon <- reconstruct_state(modSmooth$s,order_trend,order_cycle,ncycle)
interp <- interpret_state(modSmooth$s,order_trend,order_cycle,length(sc_freq),
                          flabels=sc_names)
# parms gives posterior mean estimates of parameters
parms <- sc_params_from_gibbs(sc_out)

# summarize_scout_results gives posterior statistics
# for time varying results like amplitude
df<-summarize_scout_result(sc_out,interp,out_dt=0.25,nburn=10)

# histograms of parameters
dev.new(width=6,height=9)
par(mfrow=c(4,2),mar = c(5, 5, 2, 1))
hist(sc_out$sigma2_epsilon[501:1000],xlab="",main=expression(epsilon))
hist(sc_out$sigma2_zeta[501:1000],xlab="",ylab="",main=expression(zeta))
hist(sc_out$sigma2_kappa[501:1000,1],xlab="",main=expression(kappa[1]))
hist(sc_out$sigma2_kappa[501:1000,2],xlab="",ylab="",main=expression(kappa[2]))
hist(sc_out$sigma2_kappa[501:1000,3],xlab="",main=expression(kappa[3]))
hist(sc_out$sigma2_kappa[501:1000,4],xlab="",ylab="",main=expression(kappa[4]))
hist(sc_out$sigma2_kappa[501:1000,5],xlab="Variance",main=expression(kappa[5]))
hist(sc_out$sigma2_kappa[501:1000,6],xlab="Variance",ylab="",main=expression(kappa[6]))


scparams <- sc_params_from_gibbs(sc_out)
fname <- "scha23_params.yaml"
yaml::write_yaml(scparams,fname)

write.csv(format(df,digits=3,scientific=FALSE) ,file="sc23_result_summary.csv",quote=FALSE,row.names=FALSE)


asel <- 1:18000
ts.plot(interp$amplitude$M1[asel],col='blue',ylim=c(0,1.5))
lines(interp$amplitude$M2[asel],col='black')
lines(interp$amplitude$M3[asel],col='red')
lines(interp$amplitude$M4[asel],col='purple')
lines(interp$amplitude$M5[asel],col='green')
lines(interp$amplitude$M6[asel],col='brown')
lines(interp$subtide[asel]/10)


ts.plot(sc_out$amplitude_thinned$M2[133,],ylim=c(0.,1.5))
lines(sc_out$amplitude_thinned$M1[133,],col="blue")
lines(sc_out$amplitude_thinned$M3[133,1000:16000],col="red")
lines(sc_out$amplitude_thinned$M4[133,],col="purple")
abline(h=seq(0,1.5,.5))

D1anal <- aggregate_iter(sc_out$amplitude_thinned$M1,nburn=-1)
ts.plot(D1anal,col="blue",ylim=c(0.,1.5))
D2anal <- aggregate_iter(sc_out$amplitude_thinned$M2,nburn=-1)
lines(D2anal,col="black")
D3anal <- aggregate_iter(sc_out$amplitude_thinned$M3,nburn=-1)
lines(D3anal[1000:16000],col="red")
D4anal <- aggregate_iter(sc_out$amplitude_thinned$M4,nburn=-1)
lines(D4anal,col="green")

ts.plot(interp$amplitude$M1,col='blue',ylim=c(0,1.5))
lines(interp$amplitude$M2,col='black')
lines(interp$amplitude$M3,col='red')
lines(interp$amplitude$M4,col='purple')

ts.plot(interp$subtide)
ts.plot(interp$amplitude$M6)








