# Simple test of refactoring for Freeport example

library(stochCycle)


CFS2CMS = 0.028316847
# The data are loaded in a unit that is large (cubic feet per second)
# This transformation to makes it in cubic meters and nearly unit scaled
fpt = freeport_discharge*CFS2CMS/100.

# list of tidal frequencies in radians per sample step
samples_hr = 4. # data are 15min
tide_freqs = tide_freq_radhr/samples_hr

# this subsets approximately six months,
# creates some missing data
select <- 5000:23000
y <- fpt[,1]
ysave <- fpt[select,]
#y[12300:12420] <- NA
#y[11410:11530] <- NA
#y[15000:15160] <- NA
#missing <-rbinom(n=length(y),size=1,prob=0.15)
#y[missing==1] <- NA
#ts.plot(y[1:1200])
y <- y[select]

# The D1 and D2 etc are in terms of lunar prinicpal cycle
sc_names <- c("M1","M2","M3","M4","M5","M6")
sc_freq <- tide_freqs[sc_names]
nfreq <- length(sc_freq)

mc_iter <- 1001
order_trend <- 2
order_cycle <- 3

mh_sd_rho <- 0.001 # Change this at some peril! cycle can disappear for low rho
rho_fixed <- TRUE
initial_rho <- 0.997
sigma2_zeta_fixed <- TRUE
initial_sigma2_zeta=1.e-7
initial_sigma2_kappa=1.35e-13


sc_out <- sc_gibbs(y,order_trend,order_cycle,sc_freq,mc_iter,
          rho_fixed = rho_fixed,initial_rho=initial_rho,rho_low=0.98,rho_high=0.994,mh_sd_rho=0.0005,
          initial_sigma2_epsilon=0.004,
          sigma2_zeta_fixed=sigma2_zeta_fixed,
          initial_sigma2_zeta=initial_sigma2_zeta,
          initial_sigma2_kappa)

mod <- sc_model_from_gibbs(sc_out)
modFilt <- dlm::dlmFilter(y, mod)
modSmooth <- dlm::dlmSmooth(modFilt)
recon <- reconstruct_state(modSmooth$s,order_trend,order_cycle,nfreq)
do_plot(ysave,y,recon,6000:11000)
interp <- interpret_state(modSmooth$s,order_trend,order_cycle,nfreq,flabels=c("M1","M2","M3","M4","M5","M6"))

df<-summarize_scout_result(sc_out,interp,out_dt=0.25,nburn=-1)
write.csv(format(df,digits=3,scientific=FALSE) ,file="sc23_result_summary.csv",quote=FALSE,row.names=FALSE)

# parms gives posterior mean estimates of parameters
parms <- sc_params_from_gibbs(sc_out)
yaml::write_yaml(parms,"scha_23_parms.yaml")


ts.plot(sc_out$amplitude_thinned$M2[12,],ylim=c(0.,1.5))
lines(sc_out$amplitude_thinned$M1[12,],col="blue")
lines(sc_out$amplitude_thinned$M3[12,],col="red")
lines(sc_out$amplitude_thinned$M4[12,],col="darkgreen")






