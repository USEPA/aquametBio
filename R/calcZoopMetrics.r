#' @export
#' @title Calculate zooplankton metrics for given input dataset.
#'
#' @description This function calculates all metrics associated
#' with a given input data frame using the inputs provided. These
#' inputs include the names of count, biomass, and density variables,
#' as well as a variable indicating distinctness of each taxon.
#'
#' @param inCts A data frame containing, at minimum, the variables
#' specified in the arguments for sampID, dist, ct, and taxa_id
#' and TRIBE, as well as autecology traits with names that match those
#' in the arguments ffg, habit, and ptv. In addition, there
#' should be a variable with the name in argument taxa_id that matches
#' with all of those in the indf data frame
#' @param sampID A character vector containing the names of all
#' variables in indf that specify a unique sample. If not specified,
#' the default is \emph{UID}
#' @param is_distinct A string with the name of the distinctness variable,
#' which is assumed to have only values of 0 or 1. If not specified,
#' the default is \emph{IS_DISTINCT}.
#' @param ct A string with the name of the count variable. This
#' variable should be numeric or able to be converted to numeric.
#' @param biomass A string with the name of the biomass variable. This
#' variable should be numeric or able to be converted to numeric.
#' @param density A string with the name of the density variable. This
#' variable should be numeric or able to be converted to numeric.
#' @param inTaxa a data frame containing taxonomic information,
#' including the variables PHYLUM, CLASS, SUBCLASS, ORDER, SUBORDER,
#' FAMILY, following the same taxonomy as the default taxa list,
#' \emph{zoopTaxa}. There should also be autecology traits with names
#' that match those for the arguments \emph{ffg}, \emph{clad_size},
#' and \emph{net_size}.
#' @param taxa_id A string with the name of the taxon ID variable
#' in \emph{inCts} that matches that in \emph{inTaxa}. The default
#' value is \emph{TAXA_ID}.
#' @param ffg A string with the name of the functional feeding group
#' variable in \emph{inTaxa}. Values used
#' in calculations include FILT, HERB, OMNI, PARA, PRED, representing
#' filterer, herbivore, omnivore, parasite, and predator,
#' respectively.
#' @param clad_size A string with the name of the variable in
#' \emph{inTaxa} indicating the size class of the cladoceran
#' taxon. Valid values are LARGE and SMALL.
#' @param net_size A string with the name of the variable in
#' \emph{inTaxa} indicating the net size class of a taxon.
#' Valid values are COARSE and FINE.
#' @return A data frame containing the variables in sampID and
#' the zooplankton metrics as additional variables.

calcZoopMetrics <- function(indata, sampID, is_distinct,
                           ct, biomass, density = NULL,
                           inTaxa, taxa_id='TAXA_ID',
                           ffg, clad_size, net_size){

  indata <- as.data.frame(indata)

  necVars <- c(sampID, is_distinct, ct, biomass, taxa_id)
  if(any(necVars %nin% names(indata))){
    msgTraits <- which(necVars %nin% names(indata))
    print(paste("Missing variables in input data frame:",
                paste(necVars[msgTraits], collapse=',')))
    return(NULL)
  }else{
    indata[, c(ct, biomass)] <- lapply(indata[, c(ct, biomass)], as.numeric)
  }

  if(!is.null(density)){
    if(density %nin% names(indata)){
      print("Missing density variable in input data frame.")
    }else{
      indata[, density] <- lapply(indata[, density], as.numeric)
    }
  }

  necTaxVars <- c(taxa_id, ffg, clad_size, net_size,
                  'PHYLUM', 'CLASS', 'SUBCLASS', 'ORDER',
                  'SUBORDER', 'FAMILY')
  if(any(necTaxVars %nin% names(inTaxa))){
    msgTraits <- which(necTaxVars %nin% names(inTaxa))
    print(paste("Missing variables in input taxalist:",
                paste(necTaxVars[msgTraits], collapse=',')))
    return(NULL)
  }

  # Assign characteristics to taxa
  calcData <- merge(indata, inTaxa, by = 'TAXA_ID')

  calcData$CALAN <- with(calcData, ifelse(ORDER=='CALANOIDA', 1, 0))
  calcData$COPE <- with(calcData, ifelse(SUBCLASS=='COPEPODA', 1, 0))
  calcData$COPE_HERB <- with(calcData, ifelse(COPE==1 & eval(as.name(ffg))=='HERB', 1, 0))
  calcData$ROT <- with(calcData, ifelse(PHYLUM=='ROTIFERA', 1, 0))
  calcData$COARSE <- with(calcData, ifelse(eval(as.name(net_size))=='COARSE', 1, 0))
  calcData$FINE <- with(calcData, ifelse(eval(as.name(net_size))=='FINE', 1, 0))
  calcData$SMCLAD <- with(calcData, ifelse(SUBORDER=='CLADOCERA' &
                                          eval(as.name(clad_size))=='SMALL', 1, 0))
  calcData$DAPHNIID <- with(calcData, ifelse(FAMILY=='DAPHNIIDAE', 1, 0))
  calcData$HERB <- with(calcData, ifelse(eval(as.name(ffg))=='HERB', 1, 0))
  calcData$LGCLAD <- with(calcData, ifelse(eval(as.name(clad_size))=='LARGE', 1, 0))
  calcData$OMNI <- with(calcData, ifelse(eval(as.name(ffg))=='OMNI', 1, 0))
  calcData$PLOIMA <- with(calcData, ifelse(ORDER=='PLOIMA', 1, 0))
  calcData$SIDID <- with(calcData, ifelse(FAMILY=='SIDIDAE', 1, 0))
  calcData$CRUST <- with(calcData, ifelse(PHYLUM=='ARTHROPODA' &
                          CLASS %in% c('MAXILLOPODA','BRANCHIOPODA'), 1, 0))
  calcData$CLAD <- with(calcData, ifelse(SUBORDER=='CLADOCERA', 1, 0))
  calcData$BOSM <- with(calcData, ifelse(FAMILY=='BOSMINIDAE', 1, 0))
  calcData$CYCLOP <- with(calcData, ifelse(ORDER=='CYCLOPOIDA', 1, 0))
  calcData$FLOS <- with(calcData, ifelse(ORDER=='FLOSCULARIACEAE', 1, 0))
  calcData$COLLO <- with(calcData, ifelse(ORDER=='COLLOTHECACEAE', 1, 0))
  calcData$ASPLAN <- with(calcData, ifelse(FAMILY=='ASPLANCHNIDAE', 1, 0))
  calcData$PRED <- with(calcData, ifelse(eval(as.name(ffg))=='PRED', 1, 0))
  calcData$CRUST_HERB <- with(calcData, ifelse(HERB==1 & CRUST==1, 1, 0))
  calcData$CRUST_PRED <- with(calcData, ifelse(PRED==1 & CRUST==1, 1, 0))
  calcData$CRUST_OMNI <- with(calcData, ifelse(OMNI==1 & CRUST==1, 1, 0))
  calcData$ROT_HERB <- with(calcData, ifelse(HERB==1 & ROT==1, 1, 0))
  calcData$ROT_PRED <- with(calcData, ifelse(PRED==1 & ROT==1, 1, 0))
  calcData$ROT_OMNI <- with(calcData, ifelse(OMNI==1 & ROT==1, 1, 0))
  calcData$COPE_PRED <- with(calcData, ifelse(PRED==1 & COPE==1, 1, 0))
  calcData$COPE_OMNI <- with(calcData, ifelse(OMNI==1 & COPE==1, 1, 0))
  calcData$CLAD_PRED <- with(calcData, ifelse(PRED==1 & CLAD==1, 1, 0))
  calcData$CLAD_OMNI <- with(calcData, ifelse(OMNI==1 & CLAD==1, 1, 0))
  calcData$CLAD_HERB <- with(calcData, ifelse(HERB==1 & CLAD==1, 1, 0))


}
