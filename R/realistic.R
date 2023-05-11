
# Evaluates a reduced vector to produce signal
# @param rv Reduced vector (complex)
# @param cf Central frequency
# @param t  Time
evalfreq <- function(rv,cf,t){
  out <- 0.5*(rv*exp(1i*cf*t) + Conj(rv)*exp(-1i*cf*t))
  Re(out)
}


# Create a complex admittance
cplx_admit <- function(amp,phase,Qamp,Qphase,Qexp,R0,t){
  phase <- 2*pi*phase/360.
  Qphase <- 2*pi*Qphase/360.
  R0pow <- R0^Qexp
  adjamp = amp + Qamp*R0pow
  adjphase = phase + Qphase*R0pow
  adjamp*exp(1i*adjphase)
}


#' Produce the test signal for the signal in the paper
#' @param plots Show amplitude and phase from final reduced vectors
#' @param csvfile File to write if not NULL
#' @export
realistic_tide <- function(plots=FALSE,csvfile=NULL,
                           noise_sd=0.06,noise_phi=0.12,fixed_seed=12,
                           newO1=NULL){

  set.seed(fixed_seed)

  # Construct subtide
  t = 0:(96*180)*0.25
  nt = length(t)
  noise <- colorednoise::colored_noise(nt,mean=0.,sd=noise_sd,phi=noise_phi)

  # spline describing events
  x0 = c(0.,12,27,
         34,41,43,45,46,
         47,52,58,65,
         72.,73.,74.,74.15,
         76.15,78.15,80.75,
         83,87,98,120,130,135,138,140.,
         149,150,152.,157,163,
         172,176,180)*96*0.25

  y0 = c(0.75,0.75,0.79,
         1.5,4.0,6.,7.0,6.5,
         7.5,8.5,7.,4.,
         4.7,4.,4.7,4.25,5.,
         3.,2.,3.,
         1.8,2,2.5,3.9,4.2,3.3,
         2.,3.,1.,1.2,0.7,0.4,
         1.1,1.1,0.9)
  y_0_base = spline(x0,y0,xout=t)$y
  plot(x0/24,y0)
  lines(t/24,y_0_base)

  # two low frequency cycles and an interaction
  hrs_14d = 24*14.
  hrs_180d = 24*60.
  hrs_7d = 24*7
  hrs_5d = 24*5


  d14 = 0.25*sin(2*pi*t/hrs_14d+pi)
  d180 = 0.14*sin(2*pi*t/hrs_180d+pi)
  mix = 8.*d14*d180
  d7 = 0.6*cos(2*pi*t/hrs_7d)
  d5 = 0.75*(t/10000)*cos(2*pi*t/hrs_5d)

  R0 = 0.7+(y_0_base + d14 + d180 + d7 +d5)
  R00 = 0.7+(y_0_base +d14+ d180 + d7)
  ts.plot(R0)
  lines(R00)


  d2 = cos(2*pi*t/12.5)
  d1 = cos(2*pi*t/25)


  # Diurnal and semidiurnal "reference"
  tides <- tide_freq_radhr

  amp = list()
  phase = list()

  amp[[1]] = list("2Q1"=0.02,"Q1"=0.05,"O1"=0.29,"M1"=0.01,
                  "P1"=0.015,"K1"=0.35,"J1"=0.02,"OO1"=0.012)
  phase[[1]] = list("2Q1"=95,"Q1"=97,"O1"=95,"M1"=91,
                    "P1"=102,"K1"=105,"J1"=108,"OO1"=113)


  if (! is.null(newO1)){
    amp[[1]]["O1"] <- newO1
  }


  amp[[2]] = list("2N2"=0.02,"N2"=0.16,"NU2"=0.03,"M2"=0.7,
                  "L2" =0.04,"S2"=0.18,"K2"=0.1)
  phase[[2]] = list("2N2"=290,"N2"=315.,"NU2"=321,"M2"=336,
                    "L2" =353,"S2"=336,"K2"=325)




  species=list()
  h=list()

  cunit = complex(real=0,imaginary=1.)

  C = list()
  centerfreq = c(0,0,0,0,0,0)

  print(amp[[1]])

  for (ispec in 1:2){
    C[[ispec]] = 0.
    species[[ispec]] = 0.
    flab = paste0('M',ispec)
    center = tides[[flab]]
    centerfreq[ispec] <- center
    for (fname in names(amp[[ispec]])){
      print(paste("Processing",fname))
      f=tides[[fname]]  # frequency in radians per hour
      print(f)
      print(center)
      ai = amp[[ispec]][[fname]]
      vki = f - center   # relative frequency in radians/hr
      phasei = 2*pi*phase[[ispec]][[fname]]/360.  # phase in radians
      C[[ispec]] = C[[ispec]] + ai*exp(cunit*(-phasei+vki*t))
      species[[ispec]] = species[[ispec]] + ai*cos(f*t-phasei)
      h[[ispec]] = evalfreq(C[[ispec]],center,t)
    }
  }

  for (ispec in 3:6){
    flab = paste0('M',ispec)
    center = tides[[flab]]
    centerfreq[ispec] <- center
  }



  S=list()

  sub <- 0.2*R0 + (0.2+0.2*sqrt(R0))*Mod(C[[1]]) + (0.2+0.2*sqrt(R0))*Mod(C[[2]])

  # D1 interactions

  c1a <- C[[1]]*C[[2]]*Conj(C[[2]])
  s1a <- c1a * cplx_admit(0.4,15,-0.05,-12,0.75,R0,t)

  c1b <- C[[1]]
  s1b <- c1b * cplx_admit(0.85,15,-0.1,-12,0.75,R0,t)
  S[[1]] <- s1b + s1a


  c2a <- (C[[1]]*Conj(C[[1]]))*C[[2]]
  s2a <- c2a * cplx_admit(0.1,-4,0.01,-12,0.75,R0,t)
  c2b = C[[2]]
  s2b <- c2b * cplx_admit(1.25,-5,-0.2,-12,0.75,R0,t)

  S[[2]] <-  s2a + s2b



  c3a = C[[1]]*C[[1]]*C[[1]]
  s3a = c3a* cplx_admit(0.925,-0.1,-0.12,-16,0.75,R0,t)


  c3b = C[[1]]*C[[2]]
  s3b = c3b*cplx_admit(0.925,-0.1,-0.12,-16,.75,R0,t)
  S[[3]] = s3a + s3b


  c4a = C[[2]]*C[[2]]*C[[1]]*Conj(C[[1]])
  s4a = c4a * cplx_admit(0.3,-10,-.08,-28,0.5,R0,t)

  c4b = C[[2]]*C[[2]]
  s4b = c4b * cplx_admit(0.3,-10,-.08,-28,0.5,R0,t)
  S[[4]]=s4a+s4b


  c5b = C[[2]]*C[[2]]*C[[1]]
  s5b = c5b*cplx_admit(0.08,-20,-0.008,-28,0.75,R0,t)
  S[[5]] = s5b

  #par(mfrow=c(1,1))
  #lines(Mod(s5b),col="red")

  c6a = C[[2]]*C[[2]]*C[[2]]
  s6a = c6a*cplx_admit(0.08,-20,-0.008,-35,0.5,R0,t)


  S[[6]] = s6a


  if (plots){
    par(mfrow=c(3,1))
    ts.plot(sub)
    title("Subtide")

    select <- 1:16000
    xlim = c(1,16000)

    ts.plot(Mod(S[[1]][select]),ylim=c(0,1.5),col="black")
    lines(Mod(S[[2]][select]),col="blue")
    lines(Mod(S[[3]][select]),col="red")
    lines(Mod(S[[4]][select]),col="green")
    lines(Mod(S[[5]][select]),col="purple")
    lines(Mod(S[[6]][select]),col="brown")
    title("Amplitude")


    ts.plot(Arg(S[[1]][select]),col="black",xlim=xlim,ylim=c(-4,4))
    lines(Arg(S[[2]][select]),col="blue")
    lines(Arg(S[[3]][select]),col="red")
    lines(Arg(S[[4]][select]),col="green")
    lines(Arg(S[[5]][select]),col="purple")
    lines(Arg(S[[6]][select]),col="brown")
    title("Phase")
  }

  for(ispec in 3:6){
    species[[ispec]] <- evalfreq(S[[ispec]],centerfreq[ispec],t)
  }


  signal <- sub
  for (i in 1:6){
    signal <- signal + evalfreq(S[[i]],centerfreq[i],t)
  }

  species <- do.call(cbind,species)
  colnames(species) <- c("D1","D2","D3","D4","D5","D6")

  if (!is.null(csvfile)){
    outdf <- cbind(t,signal,noise,sub,R0)
    colnames(outdf) <- c("t","signal","noise","sub","R0")
    mods <- do.call(cbind,lapply(S,Mod))
    colnames(mods) <- c("D1_amp","D2_amp","D3_amp","D4_amp","D5_amp","D6_amp")
    phases <- do.call(cbind,lapply(S,Arg))
    colnames(phases) <- c("D1_phase","D2_phase","D3_phase","D4_phase","D5_phase","D6_phase")
    outdf <- cbind(outdf,species,mods,phases)
    write.csv(format(outdf,digits=5),file=csvfile,quote=FALSE,row.names=FALSE)

  }


  out <- list(t=t,signal=signal,species=species,reduced=S,cfreq=centerfreq,sub=sub,R0=R0,noise=noise)
}
