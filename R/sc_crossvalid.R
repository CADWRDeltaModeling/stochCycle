
#' Calculate a list of blocks used for validation
#'
#' @param ylen observed data length.
#' @param foldsize size of the folds (gaps) that are deleted
#' @param blocksize spacing of the blocks. Agap is created every blocksize
#' @param first first index to include
#' @param lastmax maximum last index to include
#' @return a list of indices to drop for each step of the cross validation
#' @export
calculate_blocks <- function(ylen,foldsize=4*36,blocksize=1008,first=1,lastmax=-1)
{
  sblock <- first
  if (lastmax<0) {lastmax=ylen}
  print(nobs)
  allblocks <- list()
  iblock = 0
  while(sblock < (first+blocksize))
    {
    starts <- seq(sblock,ylen,blocksize)
    ends <- starts + foldsize - 1
    print("sblock")
    print(sblock)
    print("starts")
    print(starts)
    print(ends)
    iblock = iblock + 1
    mask <- c()
    for (istart in 1:length(starts)){
        if (ends[istart] > lastmax){next}
        newseg <- starts[istart]:ends[istart]
        mask <- c(mask,newseg)
    }
    allblocks[[iblock]] <- mask
    sblock <- sblock + foldsize
    print(paste("sblock=",sblock))
  }
   allblocks
  }

#4*36=144
#1008 is number

#' Calculate a list of blocks used for validation
#'
#' @param ylen observed data length.
#' @param foldsize size of the folds (gaps) that are deleted
#' @param blocksize spacing of the blocks. Agap is created every blocksize
#' @param first first index to include
#' @param lastmax maximum last index to include
#' @return a list of indices to drop for each step of the cross validation
#' @export
calculate_blocks2 <- function(ylen,foldsize=4*36,blocksize=1008,first=1,lastmax=-1)
{
  sblock <- first
  nobs <- length(y)
  if (lastmax<0) {lastmax=nobs}
  print(nobs)
  allblocks <- list()
  iblock = 0
  while(sblock < blocksize)
  {
    starts <- seq(sblock,ylen,blocksize)
    ends <- starts + foldsize - 1
    print("sblock")
    print(sblock)
    print("starts")
    print(starts)
    print(ends)
    iblock = iblock + 1
    mask <- c()
    for (istart in 1:length(starts)){
      if (ends[istart] > lastmax){next}
      newseg <- starts[istart]:ends[istart]
      mask <- c(mask,newseg)
    }
    allblocks[[iblock]] <- mask
    sblock <- sblock + foldsize
    print(paste("sblock=",sblock))
  }
  allblocks
}
