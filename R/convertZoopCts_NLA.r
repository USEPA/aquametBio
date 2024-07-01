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
#' @param conc_vol A string specifying the name of the numeric variable
#' with the concentrated volume of the sample in the same units as
#' tow_vol and conc_vol.
#' @param lr_taxon A string specifying the name of the variable indicating
#' a large/rare taxon. These do not have values for \emph{rawCt}.
#' This only applies to full
#' samples and not subsamples, since large/rare taxa should be ignored
#' in creating subsamples. Value should be Y or missing (NA).
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
                               tow_vol, vol_ctd, conc_vol, lr_taxon = NULL,
                               taxa_id, subsample = FALSE) {
  inCts <- as.data.frame(inCts)

  necVars <- c(
    sampID, sampType, rawCt, biofactor, taxa_id,
    tow_vol, vol_ctd, conc_vol
  )
  if (any(necVars %nin% names(inCts))) {
    msgTraits <- which(necVars %nin% names(inCts))
    print(paste(
      "Missing variables in input data frame:",
      paste(necVars[msgTraits], collapse = ",")
    ))
    return(NULL)
  }

  if (!is.null(lr_taxon)) {
    if (lr_taxon %nin% names(inCts)) {
      print("Missing variable named in lr_taxon argument - either
          change argument or add variable to input data.")

      return(NULL)
    }
  }

  if (!is.null(lr_taxon) & subsample == TRUE) {
    print("Large/rare taxa should not be included when creating
          random subsamples from the data.")
    return(NULL)
  }

  # This is includes only the non-large/rare taxa
  inCts.1 <- inCts

  inCts.1[, c(biofactor, rawCt, tow_vol, vol_ctd, conc_vol)] <-
    lapply(inCts.1[, c(biofactor, rawCt, tow_vol, vol_ctd, conc_vol)], as.numeric)
  # Drop any missing counts or any records that have any missing volume or count information
  inCts.1 <- subset(inCts.1, !is.na(eval(as.name(rawCt))))
  # |!is.na(eval(as.name(vol_ctd)))|!is.na(eval(as.name(conc_vol)))) # This does not remove those with missing biofactor or missing abundance

  inCts.1$CORR_FACTOR <- (inCts.1[, conc_vol] / inCts.1[, vol_ctd]) / inCts.1[, tow_vol]

  outCts <- inCts.1
  outCts$BIOMASS <- outCts[, rawCt] * outCts[, biofactor] * outCts$CORR_FACTOR

  outCts[which(outCts[, rawCt] == 0), "BIOMASS"] <- 0 # This indicates 0 individuals in sample, NLA taxa ID 9999
  outCts[which(is.na(outCts[, biofactor]) & outCts[, rawCt] > 0), "BIOMASS"] <- NA

  if (subsample == FALSE) {
    outCts$DENSITY <- outCts[, rawCt] * outCts$CORR_FACTOR

    outCts.1 <- aggregate(
      x = list(
        COUNT = outCts[, rawCt],
        BIOMASS = outCts$BIOMASS,
        DENSITY = outCts$DENSITY
      ),
      by = outCts[, c(sampID, sampType, taxa_id)],
      FUN = function(x) {
        sum(x, na.rm = TRUE)
      }
    )
  } else {
    outCts.1 <- aggregate(
      x = list(
        COUNT = outCts[, rawCt],
        BIOMASS = outCts$BIOMASS
      ),
      by = outCts[, c(sampID, sampType, taxa_id)],
      FUN = function(x) {
        sum(x, na.rm = TRUE)
      }
    )
  }

  # Now add back in large/rare abundances - no changes necessary, just leaves them in the data
  if (!is.null(lr_taxon) & subsample == FALSE) {
    lrCts <- subset(inCts, !is.na(eval(as.name(lr_taxon))),
      select = names(inCts) %in% c(sampID, sampType, taxa_id, lr_taxon)
    )

    lrCts.out <- lrCts

    lrCts.out <- unique(lrCts.out)

    outCts.1 <- merge(outCts.1, lrCts.out, by = c(sampID, sampType, taxa_id), all = TRUE)
  }

  outCts.1 <- outCts.1[with(outCts.1, order(
    eval(as.name(sampID)), eval(as.name(sampType)),
    eval(as.name(taxa_id))
  )), ]
  return(outCts.1)
}
