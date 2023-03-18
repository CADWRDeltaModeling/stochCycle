library(stochCycle)
library(parallel)

CFS2CMS = 0.028316847
# The data are loaded in a unit that is large (cubic feet per second)
# This transformation to makes it nearly unit scaled
fpt = freeport_discharge*CFS2CMS/100.

# list of tidal frequencies in radians per hour
samples_hr = 4. # data are 15min
tide_freqs = tide_freq_radhr/samples_hr

# this subsets approximately six months
select <- 5000:23000
#select <- 5000:18000  # for performance testing

#select <- 14000:17000 # used this reduced subset for development

y <- fpt[,1]
ysave <- fpt[select,]
nobs <- length(ysave)

# Others are M05(one cycle per two tidal days, inertial
sc_names <- c("M1","M2","M3","M4","M5","M6")
sc_freq <- tide_freqs[sc_names]
nfreq <- length(sc_freq)

mc_iter <- 32  # for debugging 20. Still converges pretty well
order_trend <- 2
order_cycle <- 4


initial_sigma2_kappa <- 1.35e-13

dframe<- data.frame(m=integer(),n=integer(),rho=double(),
                    sigma2_zeta=double(),cv_mse=double(),
                    filter_quality=double(),
                    sigma2_epsilon=double())
dframe2 <- data.frame(matrix(ncol=nfreq,nrow=0))
colnames(dframe2) <- paste("sigma2_kappa",sc_names,sep="_")
dframe <- cbind(dframe,dframe2)
save_col_names <- colnames(dframe)

cvblocks <- calculate_blocks(nobs,foldsize=4*25,blocksize=10*25*4+1,first=400,lastmax=-1)
nblock <- length(cvblocks)
print(paste("Number of cross validation blocks",nblock))

ygapfill <- ysave  # For dimension only
ygapfill[] <- NA   # This will be filled entirely with x-validation values

cv_kernel <- function(iblock,rho,initial_sigma2_zeta){
  ymask <- ysave
  if (iblock<=nblock){ # Block is used for cv. Iblock > nblock is full run
    mblock <- unlist(cvblocks[iblock])
    ymask[mblock] <- NA
    log_state_only <- TRUE
  }else{
    ymask <- ysave # only here to set trace
    log_state_only <- FALSE
  }

  sc_out <- tryCatch(
    {
      # this is the return value in the successful case
      sc_out <- sc_gibbs(ymask,order_trend,order_cycle,sc_freq,mc_iter,
                         sigma2_zeta_fixed=TRUE,
                         initial_sigma2_zeta=initial_sigma2_zeta,
                         initial_sigma2_kappa=initial_sigma2_kappa,
                         rho_fixed=TRUE,
                         initial_rho=rho,
                         thin_rate=1,
                         log_state_only=log_state_only)},
      error = function(condition){
        scout <- NA   # Unsuccessful case does not return a list
    }
  )
  if (! is.list(sc_out)){
    # The block has failed. In a series calculation, this can be used
    # to avoid going any further. In parallel, not much can be done
    success <- FALSE
    message("Failed on block")
    return(NULL)
  }
  if(iblock <= nblock){
     # block is being used for cross-validation
     # so only the filled y are needed
     nburn <- mc_iter %/% 2
     posty <- colMeans(sc_out$state_thinned[(nburn+1):mc_iter,])
     retval <- list(block=mblock,posty=posty)
  }else{
    # block is being used to estimate posterior parameters. Return everything
    retval <- list(sc_out=sc_out,block=NULL)
  }
}

in_parallel <- TRUE  # toggling off is best for debugging
if (in_parallel){
  clust <- makeCluster(min(detectCores(),length(cvblocks)))
  clusterEvalQ(clust,library("stochCycle"))

  # These are the variables that are not modified in the loop
  clusterExport(clust,c("ysave","order_trend","order_cycle",
                      "sc_freq","mc_iter",
                      "initial_sigma2_kappa","cvblocks","nblock"))
}

# Add an n+1 block that is assigned to NULL. This is a signal
# to do a full run

cvblocks <- append(cvblocks,list(NULL))

#trial_zeta <- c(5.e-8,1.e-7,5.e-7,1.e-6,5.e-6,1e-5)
third = 1./3.
trial_zeta <- 10^(seq(-7.- third,-4.9,third))
trial_rho <- c(0.991,0.992,0.993,0.994,0.995,0.996)
for(sz in trial_zeta)
{
  initial_sigma2_zeta = sz
  for (rho in trial_rho){
    print(paste("Case: rho=",rho,"sigma2_zeta=",sz))
	  success <- TRUE
    nblock_augment <- nblock+1
    if(in_parallel){
	    parallel_result <- parLapply(clust,1:nblock_augment,cv_kernel,
	                                 rho,initial_sigma2_zeta)
    }else{
      parallel_result <- lapply(1:nblock_augment,cv_kernel,
                                   rho,initial_sigma2_zeta)
    }
    for (iblock in 1:nblock){
      block <- unlist(cvblocks[iblock])
      cv_result <- parallel_result[[iblock]]
      if (is.null(cv_result)){
        success = FALSE
        break
      }
      ygapfill[block] <- cv_result$posty[cv_result$block]
    }
	  nburn <- mc_iter %/% 2
    if(success) {
      resid <- ysave-ygapfill
      mse <- mean(resid^2.,na.rm=TRUE)
      print(paste("mse",mse,"="))
    }else{
      resid<-NA
      mse <- NA
      print("mse calculation failed")
    }
    if (success){
      sc_out <- parallel_result[[nblock+1]]$sc_out
      print(sc_out$sigma2_epsilon)
      filter_eval <- regular_filter_targets_1thru4()
      filter_eval$filter_scale <- 0.004
      params <-sc_params_from_gibbs(sc_out,nburn=-1,filter_eval=filter_eval)

    }else{
      params <- list(filter_prior_value=NA,sigma2_eps=NA,sigma2_kappa=rep(NA,nfreq))
    }

    dframe <- rbind(dframe,c(order_trend,order_cycle,rho,initial_sigma2_zeta,mse,
                              params$filter_prior_val,params$sigma2_eps,
                              params$sigma2_kappa))

    colnames(dframe) <- save_col_names
    write.csv(format(dframe,digits=5),"out_sc.csv",
              quote=FALSE,row.names=FALSE)
  }
}

if (in_parallel){stopCluster(clust)}
