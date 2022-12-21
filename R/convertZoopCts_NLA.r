#' @export
#' @title Convert raw zooplankton counts to density
#' and biomass
#'
#' @description This function is only used to calculate
#' zooplankton density and biomass values from raw count
#' data using tow volume, concentrated volume, volume
#' counted, and the biomass factor. The biomass factor is
#' provided by the laboratory to convert counts to biomass.
#'
#' @param inCts A data frame containing, at minimum, the variables
#' specified in the arguments for sampID, sampType, rawct, taxa_id,
#' biofactor, tow_vol, vol_ctd, and conc_vol.
#' @param sampID A character vector containing the names of all
#' variables in indf that specify a unique sample. If not specified,
#' the default is \emph{UID}.
#' @param sampType A string specifying the name of the variable that
#' indicates the sample type.
#' @param rawCt A string specifying the name of the variable for
#' the raw count values for each taxon.
#' @param biofactor A string specifying the name of the numeric variable
#' for the biomass factor for each taxon, as provided by the laboratory.
#' @param tow_vol A string specifying the name of the numeric
#' variable with the tow volume in same units as vol_ctd and conc_vol.
#' @param vol_ctd A string specifying the name of the numeric variable
#' with the volume of the sample counted in the same units as tow_vol
#' and conc_vol.
#' @param conc_vol A strings specifying the name of the numeric variable
#' with the concentrated volume of the sample in the same units as
#' tow_vol and conc_vol.
#' @param taxa_id A string specifying the name of the variable indicating
#' the taxon being used.
#' @param subsample A logical operator to indicate whether the rawCt values
#' are part of a random subsample or the original sample. Valid values
#' are TRUE or FALSE. Default value is FALSE.
#' @return A data frame containing the \emph{sampID} and
#' \emph{sampType}, and \emph{taxa_id} fields, plus
#' BIOMASS, DENSITY, and COUNT fields.
#' @author Karen Blocksom \email{Blocksom.Karen@epa.gov}
convertZoopCts_NLA <- function(inCts, sampID, sampType, rawCt, biofactor,
                            tow_vol, vol_ctd, conc_vol, taxa_id, subsample=FALSE){
  inCts <- as.data.frame(inCts)

  necVars <- c(sampID, sampType, rawCt, biofactor, taxa_id,
               tow_vol, vol_ctd, conc_vol, subsample)
  if(any(necVars %nin% names(inCts))){
    msgTraits <- which(necVars %nin% names(inCts))
    print(paste("Missing variables in input data frame:",
                paste(necVars[msgTraits], collapse=',')))
    return(NULL)
  }

  inCts[, c(biofactor, rawCt, tow_vol, vol_ctd, conc_vol)] <-
    lapply(inCts[, c(biofactor, rawCt, tow_vol, vol_ctd, conc_vol)], as.numeric)
  # Drop any missing counts or any records that have any missing volume or biofactor information
  inCts <- subset(inCts, !is.na(eval(as.name(rawCt))) & !is.na(eval(as.name(tow_vol))) & !is.na(eval(as.name(vol_ctd))) &
                    !is.na(eval(as.name(biofactor))) & !is.na(eval(as.name(conc_vol))))

  inCts$CORR_FACTOR <- (inCts[, conc_vol]/inCts[, vol_ctd])/inCts[, tow_vol]

  outCts <- inCts
  outCts$BIOMASS <- outCts[, rawCt] * outCts[, biofactor]


  if(subsample == FALSE){
    outCts$DENSITY <- outCts[, rawCt] * outCts$CORR_FACTOR

    outCts.1 <- aggregate(x = list(COUNT = outCts[, rawCt],
                                 BIOMASS = outCts$BIOMASS,
                                 DENSITY = outCts$DENSITY),
                        by = outCts[, c(sampID, sampType, taxa_id)],
                        FUN = sum)
  }else{
    outCts.1 <- aggregate(x = list(COUNT = outCts[, rawCt],
                                   BIOMASS = outCts$BIOMASS),
                          by = outCts[, c(sampID, sampType, taxa_id)],
                          FUN = sum)
  }

  return(outCts.1)
}
