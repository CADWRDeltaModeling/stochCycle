# Simple test of refactoring for Freeport example

library(stochCycle)


sf = read.csv("noaa_sffpx_9414290_elev_2012.csv")[,2]
sf = sf[seq(1,23000/4*10+1,5)]   # every 30min



# list of tidal frequencies in radians per sample step
samples_hr = 2. # data are 15min
tide_freqs = tide_freq_radhr/samples_hr

# this subsets approximately six months,
# creates some missing data
select <- 5000:23000
select <- 2499:11501
y <-sf[select]

#y < y + 1.e-3*rnorm(length(y))
ysave <- y*1.


# The D1 and D2 etc are in terms of lunar prinicpal cycle
sc_names <- c("M1","M2","M3","M4")
sc_freq <- tide_freqs[sc_names]
nfreq <- length(sc_freq)

mc_iter <- 501
order_trend <- 2
order_cycle <- 4

rho_fixed <- TRUE
initial_rho <- 0.9935
sigma2_zeta_fixed <- FALSE
initial_sigma2_zeta=1.e-8
initial_sigma2_kappa=1.35e-13

set.seed(201) #most paper examples had seed(201)
sc_out <- sc_gibbs(y,order_trend,order_cycle,sc_freq,mc_iter,
          rho_fixed = rho_fixed,initial_rho=initial_rho,
          initial_sigma2_epsilon=0.0002,
          sigma2_zeta_fixed=sigma2_zeta_fixed,
          initial_sigma2_zeta=initial_sigma2_zeta,
          initial_sigma2_kappa)

mod <- sc_model_from_gibbs(sc_out)
modFilt <- dlm::dlmFilter(y, mod)
modSmooth <- dlm::dlmSmooth(modFilt)

recon <- reconstruct_state(modSmooth$s,order_trend,order_cycle,nfreq)
do_plot(ysave,y,recon,3000:5500)
interp <- interpret_state(modSmooth$s,order_trend,order_cycle,nfreq,flabels=c("D1","D2","D3","D4"))

asel <- 1:length(y)
ts.plot(interp$amplitude$D1,col='blue',ylim=c(0,2.5))
lines(interp$amplitude$D2,col='dark gray')

sf_d1_phase_post_mean <- aggregate_iter(sc_out$phase_thinned[['M1']],-1,is_angle=TRUE)
sf_d2_phase_post_mean <- aggregate_iter(sc_out$phase_thinned[['M2']],-1,is_angle=TRUE)

ts.plot(sc_out$amplitude_thinned$M2[31,],ylim=c(0.,.5),col='brown')
lines(sc_out$amplitude_thinned$M1[31,],col="blue")
lines(sc_out$amplitude_thinned$M3[31,],col="red")
lines(sc_out$amplitude_thinned$M4[31,],col="darkgreen")



