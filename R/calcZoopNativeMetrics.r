#' @export
#' @title Calculate percent native metrics
#' @description This function calculates percent native for
#' each of the variables in \emph{inputVars}.
#' @param indata Input data frame containing variables as identified
#' in the arguments for \emph{sampID} and \emph{inputVars}.
#' @param sampID sampID A character vector containing the names of all
#' variables in indf that specify a unique sample. If not specified,
#' the default is \emph{UID}.
#' @param inputNative A character vector with the names of
#' metrics with the totals of native taxa for various
#' characteristics (e.g., count, density, biomass, 300-count
#' subsample).
#' @param inputTotals A character vector with the names of
#' metrics with the totals of all taxa. The order must correspond
#' to the order of metrics listed in \emph{inputNative}.
#' @param nonnative A string with the name of the numeric variable
#' indicating a non-native taxon. A value of 1 indicates that a
#' taxon is non-native. All other values will be ignored.
#'
#' @return A data frame containing percent native metrics
#' for each of the \emph{inputNative} variables as a subset
#' of the \emph{inputTotals}. The names of the resulting
#' metrics are based on the names of the inputs. Metrics
#' with suffixes are replaced as follows: NIND = PIND,
#' NTAX = PTAX, DEN = PDEN, BIO = PBIO.
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
