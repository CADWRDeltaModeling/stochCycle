library(stochCycle)
require(VMDecomp)
library(pracma)
library(seewave)
require(colorednoise)





#######################
realtide <- realistic_tide(plot=TRUE)

signal <- realtide[['signal']]
#yy <- realtide[['']]
S <- realtide[['reduced']]
t = realtide[['t']]
species = realtide[['species']]

cfreq <- realtide$cfreq
sample_frq = cfreq
ny <- length(signal)
set.seed(13)
noise <- colored_noise(ny,mean=0.,sd=0.05,phi=0.1)
#noise <- 2.e-2*rnorm(length(y))
yy <- signal + noise

ts.plot(yy[1:600])


K <- 18 # res_k
res_1d <- vmd(data = yy,
              alpha = 2000,
              tau = 0,
              K = K,
              DC = FALSE,
              init = 1,
              tol = 1e-6,
              verbose = FALSE)


#spectral analysis
op <- par(mfrow = c(3,2))
sampling_interval_day <- 15.0/(24*60) # sampling duration is 15 minutes
par(mfrow=c(3,2))
for (ispec in 1:6) {
  print("item")
  item = ispec #sorder[ispec]
  print(item)
  x.spec <- spectrum(res_1d$u[,item], log="no", span=10,plot=FALSE)
  spx <- x.spec$freq/sampling_interval_day
  spy <- 2*x.spec$spec
  item_mode = glue::glue("Mode {item}")
  plot(spy~spx,xlab="frequency (cpd)",main = item_mode, ylab="spectral density",type="l", xlim=c(0,10))
}

if (K==14){
  sorder <- c(5,2,3,4)
  subindex <- c(1)
}
if (K==18){
  sorder <- c(4,2,5,3)
  subindex <- c(1)
}

indexes <- 500:1500
op <- par(mfrow = c(2,2))
for (isp in 1:4 ) {
  item = sorder[isp]
  item_mode = glue::glue("IMF {item}")
  plot(x = res_1d$u[indexes, item], type = 'l',
       main = item_mode, xlab = 'Time', ylab = '')
  lines(species[[isp]][indexes],col="red")
}


#central frequency
central_freq <- tail(res_1d$omega,1)/sampling_interval_day
central_period <- 1/central_freq
print ("central frequency (cpd)")
print (central_freq)
print ("central period (day)")
print (central_period)


subndx <- 1
envelop <- as.data.frame(cbind(t,res_1d$u[,subndx]))
names(envelop) <- c("time","subtide")

par(mfrow=c(2,2))
for (mode in 1:4){
  item <- sorder[mode]
  item_mode = glue::glue("Mode {item}")
  hilb <- hilbert(res_1d$u[,item],sample_frq)
  phase <- as.vector(angle(hilb))


  amp <- as.vector(abs(hilb))

  #phase <- rm_reference_phase(phase,refphase) + refphase
  refamp = Mod(S[[mode]])
  species_title <- paste0("D",mode)
  #ts.plot(amp)

  #lines(refamp,col="red")
  pfreq = 0.
  refphase <- cfreq[mode]*t
  compare <- Arg(S[[mode]])
  diffphase = unwrap(unwrap(phase))-refphase + pfreq*t
  diffphase[diffphase < -30] <- diffphase[diffphase < -30] + 10*pi
  diffphase[diffphase < -10] <- diffphase[diffphase < -10] + 4*pi
  if (mode %in% c(1,3)){
    compare[compare > 0] = compare[compare>0] - 2*pi
  }
  diffphase[diffphase >0]
  if (mode==1){
      ts.plot(diffphase,ylim=c(-4,0))
  }
  else{
    ts.plot(diffphase)
  }
  lines(compare,col="red")


  title(species_title)
  newcol <- cbind(amp,phase)
  newcols <- data.frame(newcol)
  names(newcols) <- paste0(species_title,c("_amp","_phase"))
  envelop <- cbind(envelop, newcols)
}

write.csv(format(envelop,digits=3,scientific=FALSE) ,file="realistic_vmd.csv",quote=FALSE,row.names=FALSE)

par(mfrow=c(1,1))
plot(t/(24),sub,type='l')
lines(t/24,envelop[,'subtide'],col='red')












