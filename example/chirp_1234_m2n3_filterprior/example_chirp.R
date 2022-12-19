# This is the SCHA(2,3) case in the paper

library(stochCycle)
set.seed(15)
order_trend <- 2
order_cycle <- 3

chirp <- as.data.frame(jay_flinchem_chirptest())
y <- chirp[,"D"]
y <- y + 2.e-2*rnorm(length(y))

# list of tidal frequencies in radians per hour
samples_hr = 4. # data are 15min
tide_freqs = tide_freq_radhr/samples_hr

# this subsets approximately six months, minimal missing data
# create some missing data
ysave <- y
select <- 1:length(y)


regress_names <- c("K1","M2","O1","S2","Q1","N2","L2","M4","MK3","MO3")
regress_freq <- tide_freqs[regress_names]

# Others are M05(one cycle per two tidal days, inertial
sc_names <- c("M1","M2","M3","M4")
sc_freq <- tide_freqs[sc_names]

# Combo for paper: 4,4,0.993, start with 1.35e-13 for both
mc_iter <- 3001
regress_names <- c("K1","M2","O1","S2","Q1","N2","L2","M4","MK3","MO3") #,"M6","MK5")

# prior for eps will be uniform
rho_initial <- 0.995
# This enables an alternate prior that is concentrated at high values
# It isn't parameterized yet
use_rho_prior <- 0.
regress <- FALSE
regress_lmr <- FALSE
sigma2_zeta_fixed <- FALSE   # previously not fixed
#1 0.001 works for flow(3,3), becomes smaller 0.1 and (2,2)
initial_sigma2_zeta <- 2.e-07   # from successful 2-4 field example
initial_sigma2_kappa <- 1.e-13 # previously 1.35e-11 for 2-4
sigma2_diffuse <- 1. #
thin_rate <- 5

initial_sigma2_epsilon <- 4.e-4
zeta_prior_shape=0.
zeta_prior_rate <- 5e-5*zeta_prior_shape*0.

# 0.005 or even 0.002 for experiment with narrow
# 0.1 for (2,4) order and rho with prior in U(0.3, .999)
# (0.15 worked for flow (3,3)
mh_sd_rho <- 0.0001
# For 2-3
mh_sd_kappa <- 1.e-14
mh_sd_zeta <- 8.e-8
mh_sd_eps <- 1.e-5
initial_sigma2_zeta <- 2.e-07   # from successful 2-4 field example
initial_sigma2_kappa <- 1.5e-13 # previously 1.35e-11 for 2-4

#all_test_freq <- regular_filter_targets_more1()
all_test_freq <- regular_filter_targets_23()
test_freqs <-  all_test_freq$test_freqs
test_targs <- all_test_freq$test_targets
filter_scale <- 0.02 # 0.05 was good for 2-4 but not for 2-3
max_rho <- 0.9985


sc_out <- sc_gibbs_filter_prior(y,order_trend,order_cycle,
          sc_freq,test_freqs,test_targs,filter_scale,samples_hr,
          regress,regress_lmr,regress_freq,initial_sigma2_epsilon,
          use_rho_prior,rho_initial,mh_sd_rho,max_rho,
          sigma2_zeta_fixed,zeta_prior_rate,zeta_prior_shape,mh_sd_zeta,mh_sd_eps,
          initial_sigma2_zeta,initial_sigma2_kappa,mh_sd_kappa,sigma2_diffuse,
          mc_iter,thin_rate)


mod <- sc_model_from_gibbs(sc_out)
modFilt <- dlm::dlmFilter(y, mod)
modSmooth <- dlm::dlmSmooth(modFilt)
recon <- reconstruct_state(modSmooth$s,order_trend,order_cycle,4)
ysave <- y
do_plot(ysave,y,recon,1:4000)
nfreq <- 4
lines(recon$full,col="red")

interp <- interpret_state(modSmooth$s,order_trend,order_cycle,length(sc_freq),
                          flabels=sc_names)

M1 <- sc_freq[["M1"]]
M2 <- sc_freq[["M2"]]
M3 <- sc_freq[["M3"]]
M4 <- sc_freq[["M4"]]

ts.plot(chirp[,"Qr"],xlim=c(3500,7500))
lines(interp$subtide,col="red")
ip1a <- aggregate_iter(sc_out$subtide_thinned)
ip1b <- aggregate_iter(sc_out$subtide_thinned,quantiles=c(0.,0.5,0.75))
lines(ip1b[1,],col="green")
lines(ip1b[2,],col="purple")
lines(ip1b[3,],col="green")
archive = list(times=chirp[,"time"],
               subtide=chirp[,"Qr"], subtide_parm=interp$subtide[2:7682],
               subtide_postmean=ip1a,
               subtide25=ip1b[1,],subtide50 = ip1b[2,],subtide75= ip1b[3,])

ip1a <- rm_reference_phase(-interp$phase$M1[2:7682],chirp[,"D1_phase"],chirp[,"D1_phasearg"])
ip1b <- rm_reference_phase(-1*sc_out$phase_thinned$M1,chirp[,"D1_phase"],chirp[,"D1_phasearg"])
ip1b <- aggregate_iter(ip1b,is_angle=TRUE)
ts.plot(ip1a,ylim=c(-1,1))
lines(chirp[,"D1_phase"],col="red")
lines(ip1b,col="green")
archive[['D1_phase']] <- chirp[,'D1_phase']
archive[['D1_phase_parm']] <- ip1a
archive[['D1_phase_postmean']] <- ip1b



ip1a <- rm_reference_phase(-interp$phase$M2[2:7682],chirp[,"D2_phase"],chirp[,"D2_phasearg"])
ip1b <- rm_reference_phase(-1*sc_out$phase_thinned$M2[400,],chirp[,"D2_phase"],chirp[,"D2_phasearg"])
ts.plot(ip1a,ylim=c(0,0.5))
lines(chirp[,"D2_phase"],col="red")
lines(ip1b,col="green")
archive[['D2_phase']] <- chirp[,'D2_phase']
archive[['D2_phase_parm']] <- ip1a
archive[['D2_phase_postmean']] <- ip1b


ip1a <- rm_reference_phase(-interp$phase$M3[2:7682],chirp[,"D3_phase"],chirp[,"D3_phasearg"])
ip1b <- rm_reference_phase(-1*sc_out$phase_thinned$M3,chirp[,"D3_phase"],chirp[,"D3_phasearg"])
ip1b <- aggregate_iter(ip1b,is_angle=TRUE)
ts.plot(ip1a,ylim=c(0,1))
lines(chirp[,"D3_phase"],col="red")
lines(ip1b,col="green")
archive[['D3_phase']] <- chirp[,'D3_phase']
archive[['D3_phase_parm']] <- ip1a
archive[['D3_phase_postmean']] <- ip1b





ip1a <- rm_reference_phase(-interp$phase$M4[2:7682],chirp[,"D4_phase"],chirp[,"D4_phasearg"])
ip1b <- rm_reference_phase(-1*sc_out$phase_thinned$M4,chirp[,"D4_phase"],chirp[,"D4_phasearg"])
ip1b <- aggregate_iter(ip1b,is_angle=TRUE)
ts.plot(ip1a,ylim=c(0,0.8))
lines(chirp[,"D4_phase"],col="red")
lines(ip1b,col="green")
archive[['D4_phase']] <- chirp[,'D4_phase']
archive[['D4_phase_parm']] <- ip1a
archive[['D4_phase_postmean']] <- ip1b





ts.plot(chirp[,"Qr"])
lines(interp$subtide,col="red")
sub <- aggregate_iter(sc_out$subtide_thinned,nburn=100,quantiles=c(0.1,0.5,0.9))
lines(sub[1,],col="green")
lines(sub[2,],col="purple")
lines(sub[3,],col="green")



ts.plot(chirp[,"D1_amp"])
lines(interp$amplitude$M1,col="red")
D1postmean <- aggregate_iter(sc_out$amplitude_thinned$M1,nburn=100)
D1anal <- aggregate_iter(sc_out$amplitude_thinned$M1,nburn=100,quantiles=c(0.1,0.5,0.9))
lines(D1anal[1,],col="green")
lines(D1anal[2,],col="purple")
lines(D1anal[3,],col="green")
archive[['D1_amp']] <- chirp[,'D1_amp']
archive[['D1_amp_parm']] <- interp$amplitude$M1[2:7682]
archive[['D1_amp_postmean']] <- D1postmean
archive[['D1_amp_25']] <- D1anal[1,]
archive[['D1_amp_50']] <- D1anal[2,]
archive[['D1_amp_75']] <- D1anal[3,]


ts.plot(chirp[,"Qr"])
lines(interp$subtide,col="red")
D2postmean <- aggregate_iter(sc_out$amplitude_thinned$M2,nburn=100)
ts.plot(chirp[,"D2_amp"])
lines(interp$amplitude$M2[2:7682],col="red")
D2anal <- aggregate_iter(sc_out$amplitude_thinned$M2,nburn=100,quantiles=c(0.1,0.5,0.9))
lines(D2anal[1,],col="green")
lines(D2anal[2,],col="purple")
lines(D2anal[3,],col="green")
archive[['D2_amp']] <- chirp[,'D2_amp']
archive[['D2_amp_parm']] <- interp$amplitude$M2[2:7682]
archive[['D2_amp_postmean']] <- D2postmean
archive[['D2_amp_25']] <- D2anal[1,]
archive[['D2_amp_50']] <- D2anal[2,]
archive[['D2_amp_75']] <- D2anal[3,]






ts.plot(chirp[,"Qr"])
lines(interp$subtide,col="red")
D3postmean <- aggregate_iter(sc_out$amplitude_thinned$M3,nburn=100)
ts.plot(chirp[,"D3_amp"],ylim=c(0.1,0.3))
lines(interp$amplitude$M3,col="red")
D3anal <- aggregate_iter(sc_out$amplitude_thinned$M3,nburn=100,quantiles=c(0.1,0.5,0.9))
lines(D3anal[1,],col="green")
lines(D3anal[2,],col="purple")
lines(D3anal[3,],col="green")
archive[['D3_amp']] <- chirp[,'D3_amp']
archive[['D3_amp_parm']] <- interp$amplitude$M3[2:7682]
archive[['D3_amp_postmean']] <- D3postmean
archive[['D3_amp_25']] <- D3anal[1,]
archive[['D3_amp_50']] <- D3anal[2,]
archive[['D3_amp_75']] <- D3anal[3,]

ts.plot(chirp[,"Qr"])
lines(interp$subtide,col="red")
D4postmean <- aggregate_iter(sc_out$amplitude_thinned$M4,nburn=100)
ts.plot(chirp[,"D4_amp"],ylim=c(0.1,0.3))
lines(interp$amplitude$M4,col="red")
D4anal <- aggregate_iter(sc_out$amplitude_thinned$M4,nburn=100,quantiles=c(0.1,0.5,0.9))
lines(D4anal[1,],col="green")
lines(D4anal[2,],col="purple")
lines(D4anal[3,],col="green")
archive[['D4_amp']] <- chirp[,'D4_amp']
archive[['D4_amp_parm']] <- interp$amplitude$M4[2:7682]
archive[['D4_amp_postmean']] <- D4postmean
archive[['D4_amp_25']] <- D4anal[1,]
archive[['D4_amp_50']] <- D4anal[2,]
archive[['D4_amp_75']] <- D4anal[3,]


ts.plot(chirp[,"Qr"])
lines(interp$subtide,col="red")
ts.plot(chirp[,"D1_amp"])
lines(interp$amplitude$M1,col="red")


ts.plot(chirp[,"D1_amp"])
lines(interp$amplitude$M1,col="red")

ts.plot(chirp[,"Qr"])
lines(interp$subtide,col="red")
ts.plot(chirp[,"D2_amp"])
lines(interp$amplitude$M2,col="red")
D2anal <- aggregate_iter(sc_out$amplitude_thinned$M2,nburn=-1)
lines(D2anal,col="green")

ts.plot(chirp[,"Qr"])
lines(interp$subtide,col="red")
ts.plot(chirp[,"D1_amp"])
lines(interp$amplitude$M1,col="red")
D1anal <- aggregate_iter(sc_out$amplitude_thinned$M1,nburn=100)
lines(D1anal,col="green")

ts.plot(chirp[,"Qr"])
lines(interp$subtide,col="red")
ts.plot(chirp[,"D3_amp"])
lines(interp$amplitude$M3,col="red")
D3anal <- aggregate_iter(sc_out$amplitude_thinned$M3,nburn=100)
lines(D3anal,col="green")

ts.plot(chirp[,"Qr"])
lines(interp$subtide,col="red")
ts.plot(chirp[,"D4_amp"])
lines(interp$amplitude$M4,col="red")
D4anal <- aggregate_iter(sc_out$amplitude_thinned$M4,nburn=-1)
lines(D4anal,col="green")


archive_df <- as.data.frame(archive)
write.csv(format(archive_df,digits=5),file="chirp23.csv",quote=FALSE,row.names=FALSE)
