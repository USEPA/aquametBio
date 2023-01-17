#' @export
#' @title Calculate percent native metrics
#' @description This function calculates percent native for
#' each of the variables in \emph{inputVars}.
#' @param indata Input data frame containing variables as identified
#' in the arguments for \emph{sampID} and \emph{inputVars}.
#' @param sampID sampID A character vector containing the names of all
#' variables in indf that specify a unique sample. If not specified,
#' the default is \emph{UID}.
#' @param is_distinct A string with the name of the distinctness variable,
#' which can have valid values of 0 (non-distinct) or 1 (distinct).
#' @param nonnative A string with the name of the numeric variable
#' indicating a non-native taxon. A value of 1 indicates that a
#' taxon is non-native. All other values will be ignored.
#'
#' @return A data frame containing percent native metrics
#' for the full data and the 300-count subsample
#'
#' @author Karen Blocksom \email{Blocksom.Karen@epa.gov}
#'
calcZoopNativeMetrics <- function(indata, sampID,
                                  inputNative, inputTotals){
  necVars <- c(sampID, inputNative, inputTotals)
  if(any(necVars %nin% names(indata))){
    msgTraits <- which(necVars %nin% names(indata))
    print(paste("Missing variables in input data frame:",
                paste(necVars[msgTraits], collapse=',')))
    return(NULL)
  }

  indata[, c(inputNative, inputTotals)] <- lapply(indata[, c(inputNative, inputTotals)], as.numeric)

  outdata <- indata

  column_names <- c(sampID, 'PARAMETER', 'RESULT')

  metsOut <- data.frame(matrix(nrow = 0, ncol = length(column_names)))
  colnames(metsOut) <- column_names

  for(i in 1:length(inputNative)){
     outdata$RESULT <- round(indata[,inputNative[i]]/indata[, inputTotals[i]]*100, 2)
     outdata$PARAMETER <- inputNative[i]

     outdata.long <- subset(outdata, select = c(sampID, 'PARAMETER', 'RESULT'))
     metsOut <- rbind(metsOut, outdata.long)
  }

  metsOut$PARAMETER <- with(metsOut, gsub('NIND', 'PIND', PARAMETER))
  metsOut$PARAMETER <- with(metsOut, gsub('NTAX', 'PTAX', PARAMETER))
  metsOut$PARAMETER <- with(metsOut, gsub('BIO', 'PBIO', PARAMETER))
  metsOut$PARAMETER <- with(metsOut, gsub('DEN', 'PDEN', PARAMETER))

  return(metsOut)
  }
