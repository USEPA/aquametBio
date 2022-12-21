#' @export
#' @title Aggregate zooplankton count data from coarse
#' and fine mesh samples into ZONW sample
#'
#' @description This function aggregates count data from NLA
#' ZOCN and ZOFN samples into a new sample type (ZONW) used in
#' metric calculations. In order to calculate the NLA
#' Zooplankton MMI, both the original and 300-count subsamples
#' for each sample type must be calculated. For the 300-count
#' subsamples, only counts and biomass apply. Repeat the function for
#' both the full samples and the 300-count subsamples.
#'
#' @param inCts A data frame containing, at minimum, the variables
#' specified in the arguments for sampID, ct, biomass, and taxa_id.
#' Include density if the dataset includes density.
#' @param inTaxa a data frame containing taxonomic information,
#' including variables for PHYLUM, CLASS, SUBCLASS, ORDER, SUBORDER,
#' FAMILY, GENUS, SPECIES, SUBSPECIES, as well as autecology traits
#' with names that match those
#' in the arguments ffg, habit, and ptv. In addition, there
#' should be a variable with the name in argument \emph{taxa_id} that matches
#' with all of those in the \emph{inCts} data frame. The default data frame
#' is zoopTaxa_nla.
#' @param sampID A character vector containing the names of all
#' variables in indf that specify a unique sample. If not specified,
#' the default is \emph{UID}.
#' @param sampType The name of the character vector containing the
#' two sample types representing fine and coarse mesh zooplankton
#' samples. Default value is \emph{SAMPLE_TYPE}.
#' @param typeFine A string containing the name of the sample type
#' for the fine mesh zooplankton sample. Default value is \emph{ZOFN}.
#' @param typeCoarse A string containing the name of the sample
#' type for the coarse mesh zooplankton sample. Default value is
#' \emph{ZOCN}.
#' @param ct A string with the name of the count variable. If not
#' specified, the default is \emph{COUNT}.
#' @param biomass A string with the name of the biomass variable.
#' If not specified, the default is \emph{BIOMASS}.
#' @param density A string with the name of the density variable.
#' If not specified, the default is \emph{NULL}.
#' @param taxa_id A string with the name of the taxon ID variable
#' in \emph{inCts} that matches that in \emph{inTaxa}. The default
#' value is \emph{TAXA_ID}.
#' @return A data frame containing the \emph{sampID} and \emph{sampType}
#' fields, plus the \emph{taxa_id},
#' \emph{biomass}, \emph{density}, and \emph{ct} fields,
#' and IS_DISTINCT. Taxonomy is aggregated to the level
#' used for the NLA MMI and distinctness is recalculated on this
#' aggregated dataset.To obtain a similar output for different
#' variations on the count, biomass, and density variables.
#' @author Karen Blocksom \email{Blocksom.Karen@epa.gov}
#'
prepZoopCombCts_NLA <- function(inCts, inTaxa=zoopTaxa_nla, sampID='UID',
                            sampType = 'SAMPLE_TYPE',
                            typeFine = 'ZOFN', typeCoarse='ZOCN',
                            ct='COUNT', biomass = 'BIOMASS',
                            density = NULL, taxa_id='TAXA_ID'){

  inCts <- as.data.frame(inCts)
  inTaxa <- as.data.frame(inTaxa)

  necVars <- c(sampID, ct, biomass, sampType, taxa_id)
  if(any(necVars %nin% names(inCts))){
    msgTraits <- which(necVars %nin% names(inCts))
    print(paste("Missing variables in input data frame:",
                paste(necVars[msgTraits], collapse=',')))
    return(NULL)
  }

  if(!is.null(density)){
    if(density %nin% names(inCts)){
      print("Missing density variable in input data frame.")
    }
  }

  necTypes <- c(typeCoarse, typeFine)
  if(any(necTypes %nin% inCts[, sampType])){
    msgType <- which(necTypes %nin% inCts[, sampType])
    print(paste("Missing sample type(s) in input data frame:",
                paste(necTypes[msgType], collapse=', ')))
  }

  inCts <- subset(inCts, select=names(inCts) %in% c(sampID, sampType, ct,
                                                    biomass, density, taxa_id))

  if(!is.null(density)){
    inCts[, c(biomass, density, ct)] <- lapply(inCts[, c(biomass, density, ct)], as.numeric)
  }else{
    inCts[, c(biomass, ct)] <- lapply(inCts[, c(biomass, ct)], as.numeric)
  }

  outCts <- aggregate(x = list(COUNT = inCts$COUNT,
                               BIOMASS = inCts$BIOMASS,
                               DENSITY = inCts$DENSITY),
                      by = inCts[, c(sampID, taxa_id)],
                      FUN = sum)

  outCts[, sampType] <- 'ZONW'

  return(outCts)

}
