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
sc_names <- c("M1","M2","M3","M4","M5","M6") #todo: added M5, was never M8
sc_freq <- tide_freqs[sc_names]

mc_iter <- 1001
regress_names <- c("K1","M2","O1","S2","Q1","N2","L2","M4","MK3","MO3") #,"M6","MK5")
order_trend <- 2
order_cycle <- 4

# Leftover
zeta_prior_shape=0.
zeta_prior_rate <- 0.

sigma2_diffuse <- 100. #
thin_rate <- 5



# This enables an alternate prior that is concentrated at high values
# It isn't parameterized yet
use_rho_prior <- 0.
regress <- FALSE
regress_lmr <- FALSE
sigma2_zeta_fixed <- FALSE   # previously not fixed
#1 0.001 works for flow(3,3), becomes smaller 0.1 and (2,2)
#initial_sigma2_zeta <- 2.5e-07   # from successful 2-4 field example
#initial_sigma2_kappa <- 1.e-13 # previously 1.35e-11 for 2-4
#initial_sigma2_epsilon <- 3.78e-4 # was 1e-3 then 4e-4


# Starting values
mh_sd_rho <- 0.0005   #0.001
mh_sd_kappa <- 5.e-14  # was 1e-14
mh_sd_zeta <- 2.e-7
mh_sd_epsilon <- 0.00002   #was 0.001
initial_rho <- 0.994
initial_sigma2_zeta <- 1.e-05   # started 1e-7
initial_sigma2_kappa <- 4.35e-13 # started 1.35e-13

max_rho <- 0.994  # bumped down


# Restart with tiny bump down in rho
initial_rho <- 0.993
max_rho <- 0.9934
initial_sigma2_epsilon <- 0.000333401480367305
initial_sigma2_zeta <- 1.02009854152571e-05
if (TRUE){
initial_sigma2_kappa <- c(1.83470322164004e-13,3.6740626572885e-12,
                          7.36955368686414e-12,2.03651869940484e-11,
                          1.3829292240856e-13,1.36567841259304e-13)
mh_sd_rho <- 0.0002   #0.001
mh_sd_kappa <- c(1e-14,1e-13,1e-13,5e-12,1e-14,1e-14)
mh_sd_zeta <- 2.e-7
mh_sd_epsilon <- 0.000005   #was 0.001
}


#all_test_freq <- regular_filter_targets_24a()
all_test_freq <- regular_filter_targets_1thru6()
test_freqs <-  all_test_freq$test_freqs
test_targs <- all_test_freq$test_targets
filter_scale <- 0.004


sc_out <- sc_gibbs_filter_prior(y,order_trend,order_cycle,
          sc_freq,test_freqs,test_targs,filter_scale,samples_hr,
          regress,regress_lmr,regress_freq,initial_sigma2_epsilon,
          use_rho_prior,initial_rho,mh_sd_rho,max_rho,
          sigma2_zeta_fixed,zeta_prior_rate,zeta_prior_shape,mh_sd_zeta,mh_sd_epsilon,
          initial_sigma2_zeta,initial_sigma2_kappa,mh_sd_kappa,sigma2_diffuse,
          mc_iter,thin_rate)

mod <- sc_model_from_gibbs(sc_out)
modFilt <- dlm::dlmFilter(y, mod)
modSmooth <- dlm::dlmSmooth(modFilt)


# Collect output
ncycle <- length(sc_freq)
interp <- interpret_state(modSmooth$s,order_trend,order_cycle,length(sc_freq),
                          flabels=sc_names)
# parms gives posterior mean estimates of parameters
parms <- sc_params_from_gibbs(sc_out)
yaml::write_yaml(parms,"scha_24_parms.yaml")
write.csv(format(df,digits=3,scientific=FALSE) ,file="sc24_result_summary.csv",quote=FALSE,row.names=FALSE)



# summarize_scout_results gives posterior statistics
# for time varying results like amplitude
df<-summarize_scout_result(sc_out,interp,out_dt=0.25,nburn=10)

# histograms of parameters
dev.new(width=6,height=9)
par(mfrow=c(4,2),mar = c(5, 5, 2, 1))
hist(sc_out$sigma2_epsilon[401:1000],xlab="",main=expression(epsilon))
hist(sc_out$sigma2_zeta[401:1000],xlab="",ylab="",main=expression(zeta))
hist(sc_out$sigma2_kappa[401:1000,1],xlab="",main=expression(kappa[1]))
hist(sc_out$sigma2_kappa[401:1000,2],xlab="",ylab="",main=expression(kappa[2]))
hist(sc_out$sigma2_kappa[401:1000,3],xlab="",main=expression(kappa[3]))
hist(sc_out$sigma2_kappa[401:1000,4],xlab="",ylab="",main=expression(kappa[4]))
hist(sc_out$sigma2_kappa[401:1000,5],xlab="Variance",main=expression(kappa[5]))
hist(sc_out$sigma2_kappa[401:1000,6],xlab="Variance",ylab="",main=expression(kappa[6]))

# Spot checking below
asel <- 1:18000
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


abline(h=c(0))




ts.plot(sc_out$amplitude_thinned$M2[33,asel],ylim=c(0.,1.5))
lines(sc_out$amplitude_thinned$M1[33,asel],col="blue")
lines(sc_out$amplitude_thinned$M3[33,asel],col="red")
lines(sc_out$amplitude_thinned$M4[33,asel],col="purple")

D1anal <- aggregate_iter(sc_out$amplitude_thinned$M1,nburn=-1)
ts.plot(D1anal,col="blue",ylim=c(0.,1.5))
D2anal <- aggregate_iter(sc_out$amplitude_thinned$M2,nburn=-1)
lines(D2anal,col="black")
D3anal <- aggregate_iter(sc_out$amplitude_thinned$M3,nburn=-1)
lines(D3anal,col="red")
D4anal <- aggregate_iter(sc_out$amplitude_thinned$M4,nburn=-1)
lines(D4anal,col="green")


ts.plot(interp$subtide)















