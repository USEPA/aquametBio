#' @export
#' @title Draw a random subsample of a specified size
#'
#' @description This function allows the user to create a random
#' subsample of a specified size from the original count data.
#' The data should be in a long format, with each taxon on a
#' separate row. The user can make this process repeatable by
#' specifying a random seed that can be used to reproduce the
#' results, as long as all of the data are in the same order.
#' This function was adapted from a script originally
#' created by John Van Sickle (rarify.r, 10 June 2005).
#'
#' @param inCts A data frame containing, at minimum, the variables
#' specified in the arguments for sampID and abund, and some
#' variable to indicate taxon, either by an ID number or name.
#' @param sampID A character vector containing the names of all
#' variables in indf that specify a unique sample. If not specified,
#' the default is \emph{UID}.
#' @param abund A string with the name of the count variable. If not
#' specified, the default is \emph{COUNT}. Any rows with a missing value
#' for the \emph{abund} variable will be dropped in the output
#' dataset.
#' @param subsize A numeric value to indicate the size of the random
#' subsample to be created.
#' @param seed Optional number to set the random seed. This allows the
#' process to be repeated and receive the same results, as long as
#' the input data are ordered the same way.
#' @return A data frame containing the same fields at the input
#' dataset, with the \emph{abund} containing the new subsample
#' count for each taxon.

#' @author Karen Blocksom \email{Blocksom.Karen@epa.gov}
#'
rarifyCounts <- function(inCts, sampID='UID', abund='COUNT',
                         subsize = NULL,
                         seed = NULL){
  start.time=proc.time()
  inCts <- subset(inCts, !is.na(eval(as.name(abund))) & as.numeric(eval(as.name(abund)))>0)

  for(i in 1:length(sampID)){
    if(i==1){
      inCts$sampid <- inCts[, sampID[i]]
    }else{
      inCts$sampid <- paste(inCts[, 'sampid'], inCts[, sampID[i]], sep='.')
    }
  }

  outCts <- inCts
  samps <- unique(outCts$sampid)
  nsamp <- length(unique(outCts$sampid))
  # parameters are set up
  # zero out all abundances in output data set
  outCts[, abund] <- 0
  # loop over samples, rarify each one in turn

  for(i in 1:nsamp) {
    # extract current sample
    isamp <- samps[i]
    flush.console()
    print(as.character(isamp))
    onesamp <- inCts[inCts$sampid==isamp,]
    onesamp <- data.frame(onesamp, row.id = seq(1, dim(onesamp)[[1]])) # add sequence numbers as a new column
    # expand the sample into a vector of individuals
    samp.expand <- rep(x=onesamp$row.id, times=onesamp[, abund])
    nbug <- length(samp.expand) # number of bugs in sample

    # If seed supplied, set random seed
    if(!is.null(seed)){
      set.seed(seed)
    }
    # vector of uniform random numbers
    ranvec <- runif(n = nbug)
    # sort the expanded sample randomly
    samp.ex2 <- samp.expand[order(ranvec)]
    # keep only the first piece of ranvec, of the desired fixed count size
    # if there are fewer bugs than the fixed count size, keep them all
    if(nbug > subsize){
      subsamp <- samp.ex2[1:subsize]
    }else{
      subsamp <- samp.ex2
      }
    # tabulate bugs in subsample
    subcnt <- table(subsamp)
    # define new subsample frame and fill it with new reduced counts
    newsamp <- onesamp
    newsamp[, abund] <- 0
    newsamp[match(newsamp$row.id, names(subcnt), nomatch = 0)>0, abund] <- as.vector(subcnt)
    outCts[outCts$sampid==isamp, abund] <- newsamp[, abund]
  }; # end of sample loop

  elaps <- proc.time() - start.time
  print(c("Rarefaction complete. Number of samples = ", nsamp), quote=F)
  print(c("Execution time (sec)= ", elaps[1]), quote=F)

  outCts$sampid <- NULL
  return(outCts)
}
