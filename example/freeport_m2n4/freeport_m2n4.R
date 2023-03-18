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

mc_iter <- 1001
order_trend <- 2
order_cycle <- 4

mh_sd_rho <- 0.001 # Change this at some peril! cycle can disappear for low rho
rho_fixed <- TRUE
initial_rho <- 0.9935
sigma2_zeta_fixed <- FALSE
initial_sigma2_zeta=1.6e-6
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
do_plot(ysave,y,recon,6000:11000)
interp <- interpret_state(modSmooth$s,order_trend,order_cycle,4,flabels=c("D1","D2","D3","D4","D5","D6"))

ts.plot(sc_out$amplitude_thinned$M2[12,],ylim=c(0.,1.5))
lines(sc_out$amplitude_thinned$M1[12,],col="blue")
lines(sc_out$amplitude_thinned$M3[12,],col="red")
lines(sc_out$amplitude_thinned$M4[12,],col="darkgreen")

# Collect output
ncycle <- length(sc_freq)
interp <- interpret_state(modSmooth$s,order_trend,order_cycle,length(sc_freq),
                          flabels=sc_names)
# parms gives posterior mean estimates of parameters
parms <- sc_params_from_gibbs(sc_out)
yaml::write_yaml(parms,"scha_24_parms.yaml")




# summarize_scout_results gives posterior statistics
# for time varying results like amplitude
df<-summarize_scout_result(sc_out,interp,out_dt=0.25,nburn=-1)
write.csv(format(df,digits=3,scientific=FALSE) ,file="sc24_result_summary.csv",quote=FALSE,row.names=FALSE)


# histograms of parameters
dev.new(width=6,height=9)
par(mfrow=c(4,2),mar = c(5, 5, 2, 1))
histselect <- 400:1000
hist(sc_out$sigma2_epsilon[histselect],xlab="",main=expression(epsilon))
hist(sc_out$sigma2_zeta[histselect],xlab="",ylab="",main=expression(zeta))
hist(sc_out$sigma2_kappa[histselect,1],xlab="",main=expression(kappa[1]))
hist(sc_out$sigma2_kappa[histselect,2],xlab="",ylab="",main=expression(kappa[2]))
hist(sc_out$sigma2_kappa[histselect,3],xlab="",main=expression(kappa[3]))
hist(sc_out$sigma2_kappa[histselect,4],xlab="",ylab="",main=expression(kappa[4]))
hist(sc_out$sigma2_kappa[histselect,5],xlab="Variance",main=expression(kappa[5]))
hist(sc_out$sigma2_kappa[histselect,6],xlab="Variance",ylab="",main=expression(kappa[6]))

# Spot checking below
asel <- 1:18000
asel <-1:4500
ts.plot(interp$amplitude$M1[asel],col='blue',ylim=c(0,1.5))
lines(interp$amplitude$M2[asel],col='dark gray')
lines(interp$amplitude$M3[asel],col='red')
lines(interp$amplitude$M4[asel],col='purple')
lines(interp$amplitude$M5[asel],col='green')
lines(interp$amplitude$M6[asel],col='brown')
lines(interp$subtide[asel]/10)

bsel <- 13000:15000
ts.plot(interp$freq$M1[bsel],col='black',ylim=c(-2.75,2.75))
lines(interp$freq$M2[bsel],col='dark grey')
lines(interp$freq$M3[bsel],col='green')



ts.plot(-interp$phase$M1[bsel],col='black',ylim=c(-3.25,3.25))
lines(-interp$phase$M2[bsel],col='dark grey')



pot2 <- read.csv("../sf_potential/sfpot_d2_post_mean.csv",header=FALSE)[,2]
pot1 <- read.csv("../sf_potential/sfpot_d1_post_mean.csv",header=FALSE)[,2]
ts.plot(pot1)
diff <- rm_reference_phase(interp$phase$D1,pot1)
