#' @export
#' @title Calculate zooplankton dominance metrics
#' @description This function calculates all of the dominance metrics used
#' for zooplankton.
#' @param df Input data frame, containing SAMPID as variable identifying unique samples
#' @param sampID A character vector containing the names of all
#' variables in indf that specify a unique sample. If not specified,
#' the default is \emph{UID}.
#' @param is_distinct A string with the name of the distinctness variable,
#' which is assumed to have only values of 0 or 1. If not specified,
#' the default is \emph{IS_DISTINCT}.
#' @param valsIn A vector with the names of numeric variables to be used
#' in calculating dominance.
#' @param valsOut A vector of output suffixes to append to dominance
#' metric, matched with the order of valsIn variables.
#' @param taxa_id A string with the name of the variable that distinctly
#' identifies taxa in each sample.
#' @param subgrp A string with the name of a subgroup in \emph{indata}, with
#' a valid value of 1 to indicate inclusion in the subgroup.
#' @return A data frame with \emph{sampID} variables and the metric containing
#' the % individuals in the
#' dominant (topN) taxa
calcZoopDomMetrics <- function(indata, sampID, is_distinct,
                               valsIn, valsOut, taxa_id,
                               subgrp = NULL){
  # convert input data to data frame
  indata <- as.data.frame(indata)
  # Check for necessary variables in dataset, if all there, convert
  # valsIn and is_distinct variables to numeric
  necVars <- c(sampID, is_distinct, valsIn, subgrp)

  if(any(necVars %nin% names(indata))){
    msgTraits <- which(necVars %nin% names(indata))
    print(paste("Missing variables in input data frame:",
                paste(necVars[msgTraits], collapse=',')))
    return(NULL)
  }else{
    indata[, c(valsIn, is_distinct)] <- lapply(indata[, c(valsIn, is_distinct)], as.numeric)
  }
  # Create empty data frame
  column_names <- c(sampID, 'PARAMETER', 'RESULT')

  outDom <- data.frame(matrix(nrow = 0, ncol = length(column_names)))
  colnames(outDom) <- column_names

  # Subset data to subgroup if needed and to only distinct values
  for(i in 1:length(valsIn)){
    if(!is.null(subgrp)){
      indata.1 <- subset(indata, eval(as.name(subgrp))==1 &
                         !is.na(eval(as.name(valsIn[i]))) &
                           eval(as.name(valsIn[i]))>0 &
                         eval(as.name(is_distinct))==1)
    }else{
      indata.1 <- subset(indata,!is.na(eval(as.name(valsIn[i]))) &
                            eval(as.name(valsIn[i]))>0 &
                            eval(as.name(is_distinct))==1)
    }
  # Call the zoopDominance function
    for(j in seq(from=1, to=5, by=2)){
      dd <- zoopDominance(indata.1, sampID, topN=j, varIn=valsIn[i], taxa_id)
      domtype <- valsOut[i]

      ee <- reshape(dd, idvar = sampID, direction = 'long',
                    varying = 'dompind', timevar = 'PARAMETER',
                    v.names = 'RESULT', times = 'dompind')

      if(is.null(subgrp)){
        ee$PARAMETER <- paste0('DOM', j, '_', domtype)
      }else{
        ee$PARAMETER <- paste0('DOM', j, '_', subgrp, '_', domtype)
      }

      outDom <- rbind(outDom, ee)
    }
  }
  return(outDom)
}
