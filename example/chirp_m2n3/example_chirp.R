# This is the SCHA(2,4) case in the paper

library(stochCycle)
set.seed(13)
chirp <- as.data.frame(jay_flinchem_chirptest())
y <- chirp[,"D"]
y <- y + 2.e-2*rnorm(length(y))

# list of tidal frequencies in radians per hour
samples_hr = 4. # data are 15min
tide_freqs = tide_freq_radhr/samples_hr

sc_names <- c("M1","M2","M3","M4")
sc_freq <- tide_freqs[sc_names]



mc_iter <- 3001
order_trend <- 2
order_cycle <- 3

rho_fixed <- TRUE
initial_rho <- 0.997
 sigma2_zeta_fixed <- FALSE
initial_sigma2_zeta=1.e-5 #6.5e-6
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
interp <- interpret_state(modSmooth$s,order_trend,order_cycle,length(sc_freq),
                          flabels=sc_names)

recon <- reconstruct_state(modSmooth$s,order_trend,order_cycle,4)
ysave <- y
do_plot(ysave,y,recon,1:4000)
nfreq <- 4
lines(recon$full,col="red")




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
D1anal <- aggregate_iter(sc_out$amplitude_thinned$M1,nburn=-1)
lines(D1anal,col="green")

ts.plot(chirp[,"Qr"])
lines(interp$subtide,col="red")
ts.plot(chirp[,"D3_amp"],ylim=c(.1,.3))
lines(interp$amplitude$M3,col="red")
D3anal <- aggregate_iter(sc_out$amplitude_thinned$M3,nburn=-1)
lines(D3anal,col="green")

ts.plot(chirp[,"Qr"])
lines(interp$subtide,col="red")
ts.plot(chirp[,"D4_amp"],ylim=c(0,.15))
lines(interp$amplitude$M4,col="red")
D4anal <- aggregate_iter(sc_out$amplitude_thinned$M4,nburn=-1)
lines(D4anal,col="green")


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

ip1a <- rm_reference_phase(-interp$phase$M1[2:7682],
                           chirp[,"D1_phasearg"])+chirp[,"D1_phase"]
ip1b <- rm_reference_phase(-1*sc_out$phase_thinned$M1,
                           chirp[,"D1_phasearg"])
ip1b <- sweep(ip1b,2,chirp[,"D1_phase"],"+")
ip1b <- aggregate_iter(ip1b,is_angle=TRUE)
ts.plot(ip1a,ylim=c(-0.2,0.4))
lines(chirp[,"D1_phase"],col="red")
lines(ip1b,col="green")
archive[['D1_phase']] <- chirp[,'D1_phase']
archive[['D1_phase_parm']] <- ip1a
archive[['D1_phase_postmean']] <- ip1b



ip1a <- rm_reference_phase(-interp$phase$M2[2:7682],
                           chirp[,"D2_phasearg"])+chirp[,"D2_phase"]
ip1b <- rm_reference_phase(-1*sc_out$phase_thinned$M2,
                           chirp[,"D2_phasearg"])
ip1b <- sweep(ip1b,2,chirp[,"D2_phase"],"+")
ip1b <- aggregate_iter(ip1b,is_angle=TRUE)
ts.plot(ip1a,ylim=c(0,0.5))
lines(chirp[,"D2_phase"],col="red")
lines(ip1b,col="green")
archive[['D2_phase']] <- chirp[,'D2_phase']
archive[['D2_phase_parm']] <- ip1a
archive[['D2_phase_postmean']] <- ip1b


ip1a <- rm_reference_phase(-interp$phase$M3[2:7682],
                           chirp[,"D3_phasearg"]) +chirp[,"D3_phase"]
ip1b <- rm_reference_phase(-1*sc_out$phase_thinned$M3,
                           chirp[,"D3_phasearg"])
ip1b <- sweep(ip1b,2,chirp[,"D3_phase"],"+")
ip1b <- aggregate_iter(ip1b,is_angle=TRUE)
ts.plot(ip1a,ylim=c(0,1))
lines(chirp[,"D3_phase"],col="red")
lines(ip1b,col="green")
archive[['D3_phase']] <- chirp[,'D3_phase']
archive[['D3_phase_parm']] <- ip1a
archive[['D3_phase_postmean']] <- ip1b





ip1a <- rm_reference_phase(-interp$phase$M4[2:7682],
                           chirp[,"D4_phasearg"]) + chirp[,"D4_phase"]
ip1b <- rm_reference_phase(-1*sc_out$phase_thinned$M4,
                           chirp[,"D4_phasearg"])
ip1b <- sweep(ip1b,2,chirp[,"D3_phase"],"+")
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
ts.plot(chirp[,"D4_amp"],ylim=c(0.05,0.13))
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

