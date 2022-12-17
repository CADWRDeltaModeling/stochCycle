


#' Jay(1999) tide with subtide interference
#' @param n cycles per day in subtide
#' @export
jay_subtide_interfere <- function(n=4)
{
  samples_hr = 4. # data are 15min
  tide_freqs = tide_freq_radhr/samples_hr
  print(tide_freqs)
  M1 <- tide_freqs["M1"]
  M2 <- tide_freqs["M2"]
  M3 <- tide_freqs["M3"]
  t <- 1:(805*samples_hr)
  cos(M1*t)+sin(M1*t) +2.5*cos(M2*t)+2.5*sin(M2*t)+cos(M3*t)+sin(M3*t)+
      4.*cos(2.*3.14159*t/(24*samples_hr*n))

}


#' Jay(1999) tide with subtide interference
#'
#' @param n cycles per day in subtide
#' @param freqs vector of frequencies to include (M2,O1,K1,SUB)
#' @export
jay2 <- function(n=4, freqs=c("M2","O1","K1","SUB"))
{
  samples_hr = 4. # data are 15min
  tide_freqs = tide_freq_radhr/samples_hr
  print(tide_freqs)
  K1 <- tide_freqs["K1"]
  O1 <- tide_freqs["O1"]
  M2 <- tide_freqs["M2"]
  M3 <- tide_freqs["M3"]
  f <- tide_freqs["f"]

  t <- 1:(805*samples_hr)
  res <- 0.*t
  if("SUB" %in% freqs){
    res <- res + 4.*cos(2.*pi*t/(24*samples_hr*n))
  }
  if("M2" %in% freqs){
    res <- res + 0.575*cos(M2*t-5.9)
  }
  if("O1" %in% freqs){
    res <- res+0.23*cos(O1*t-1.89)
  }
  if("K1" %in% freqs){
    res <- res +0.37*cos(K1*t-1.83)
  }
  if("D1D2" %in% freqs){
    res <- res + (0.23*cos(O1*t-1.89)*0.37*cos(K1*t-1.83))
  }
  if("M3" %in% freqs){
    res <- res +0.25*cos(M3*t)
  }
  if("f" %in% freqs){
    res <- res +1.0*cos(f*t)
  }
  res
}


#' Jay(1999) tide with subtide interference
#'
#' @param n cycles per day in subtide
#' @param freqs vector of frequencies to include (M2,O1,K1,SUB)
#' @export
jay3 <- function(n=4, freqs=c("M2","O1","K1","SUB"),nhour=805,samples_hr=4.)
{
  tide_freqs = tide_freq_radhr/samples_hr
  print(tide_freqs)
  K1 <- tide_freqs["K1"]
  O1 <- tide_freqs["O1"]
  M1 <- tide_freqs["M1"]
  M2 <- tide_freqs["M2"]
  M3 <- tide_freqs["M3"]
  f <- tide_freqs["f"]

  t <- 1:(nhour*samples_hr)
  res <- 0.*t
  #if("SUB" %in% freqs){
  #  res <- res + 4.*cos(2.*pi*t/(24*samples_hr*n))
  #}

  if("M2" %in% freqs){
    res <- res + 1.*cos(M2*t)
  }
  if("O1" %in% freqs){
    res <- res+0.35*cos(O1*t)
  }
  if("K1" %in% freqs){
    res <- res +0.35*cos(K1*t)
  }
  if("f" %in% freqs){
    res <- res +1.0*cos(f*t)
  }
  res
}



#modulate_inertial(f){
#  f
#}

