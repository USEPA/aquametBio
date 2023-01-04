#' @export
#' @title Calculate zooplankton totals
#' @description This internal function calculates totals for
#' counts, biomass, and density from inputs. This function is
#' called from within various metric calculation functions using
#' various data subsets.
#' @param indata Input data frame containing variables as identified
#' in the arguments for \emph{sampID} and \emph{inputNames}.
#' @param sampID sampID A character vector containing the names of all
#' variables in indf that specify a unique sample. If not specified,
#' the default is \emph{UID}.
#' @param is_distinct A string with the name of the distinctness variable,
#' which can have valid values of 0 (non-distinct) or 1 (distinct).
#' @param inputSums A vector with the names of numeric variables to be
#' summed. If variables are not numeric, they will be converted.
#' This may lead to NAs being produced if any values cannot be
#' converted to numeric.
#' @param outputSums A vector with the names of numeric variables
#' containing total of each input variable.
#' @param outputTaxa A string with the name of the variable for the
#' number of taxa output variable.
#'
#' @return A data frame containing the totals of input variables by
#' the variables in \emph{sampID}.
#'
#' @author Karen Blocksom \email{Blocksom.Karen@epa.gov}
#'
calcZoopTotals <- function(indata, sampID, is_distinct,
                           inputSums, outputSums,
                           outputTaxa){
  necVars <- c(sampID, inputSums)
  if(any(necVars %nin% names(indata))){
    msgTraits <- which(necVars %nin% names(indata))
    print(paste("Missing variables in input data frame:",
                paste(necVars[msgTraits], collapse=',')))
    return(NULL)
  }

  indata[, c(inputSums, is_distinct)] <- lapply(indata[, c(inputSums, is_distinct)], as.numeric)

  outdata <- unique(indata[, sampID])

  for(i in 1:length(inputSums)){

    temp <- aggregate(x=list(x.out = indata[, inputSums[i]]),
                                   by = indata[, c(sampID)],
                                   FUN = function(x){sum(x, na.rm=TRUE)})

    outdata <- merge(outdata, temp)
    names(outdata)[names(outdata)==inputSums[i]] <- outputSums[i]
  }

  # now count distinct taxa using distinctness variable
  ntax <- aggregate(x = list(ntax <- indata[, is_distinct]),
                    by = indata[, sampID],
                    FUN = function(x){sum(x, na.rm=TRUE)})

  names(ntax)[names(ntax) == is_distinct] <- outputTaxa

  outdata <- merge(outdata, ntax, by = sampID)

  return(outdata)
}
