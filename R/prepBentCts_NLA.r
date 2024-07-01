#' @export
#' @title Aggregate benthic macroinvertebrate data for use
#' with NLA benthic metrics and MMIs
#'
#' @description This function aggregates count data to match
#' those taxonomic levels used in NLA and in the NLA MMI
#'
#' @param inCts A data frame containing, at minimum, the variables
#' specified in the arguments for sampID, ct, and taxa_id.
#' @param inTaxa a data frame containing taxonomic information,
#' including variables for PHYLUM, CLASS, ORDER, FAMILY, SUBFAMILY,
#' and TRIBE, as well as autecology traits with names that match those
#' in the arguments ffg, habit, and ptv. In addition, there
#' should be a variable with the name in argument taxa_id that matches
#' with all of those in the indf data frame. The default data frame
#' is bentTaxa_nla.
#' @param sampID A character vector containing the names of all
#' variables in indf that specify a unique sample. If not specified,
#' the default is \emph{UID}.
#' @param ct A string with the name of the count variable. If not
#' specified, the default is \emph{TOTAL}.
#' @param taxa_id A string with the name of the taxon ID variable
#' in \emph{inCts} that matches that in \emph{inTaxa}. The default
#' value is \emph{TAXA_ID}.
#' @return A data frame containing the \emph{sampID} fields, plus
#' TAXA_ID, TOTAL, and IS_DISTINCT. Taxonomy is aggregated to the level
#' used for the NLA MMI and distinctness is recalculated on this
#' aggregated dataset.
#' @author Karen Blocksom \email{Blocksom.Karen@epa.gov}
#'
prepBentCts_NLA <- function(
    inCts, inTaxa = bentTaxa_nla, sampID = "UID", ct = "TOTAL",
    taxa_id = "TAXA_ID") {
  # Convert data into data frames just in case
  inCts <- as.data.frame(inCts)
  inTaxa <- as.data.frame(inTaxa)

  ctVars <- c(sampID, ct, taxa_id)
  if (any(ctVars %nin% names(inCts))) {
    msgTraits <- which(ctVars %nin% names(inCts))
    print(paste("Missing variables in input data frame:", paste(ctVars[msgTraits], collapse = ",")))
    return(NULL)
  }

  inCts <- subset(inCts, select = c(sampID, ct, taxa_id))

  # If necessary, load the bentTaxa data frame and assign it to inTaxa.  Though
  # NON_TARGET taxa are included in the table provided by NRSA, we need to exclude
  # them from our calculations.
  if (is.null(inTaxa)) {
    inTaxa <- bentTaxa_nla
    inTaxa <- subset(inTaxa, is.na(NON_TARGET) | NON_TARGET == "")
  }

  ## This code assumes that the following are columns in the taxa file: PHYLUM, CLASS, ORDER, FAMILY, SUBFAMILY, TRIBE, HABIT, FFG, PTV,
  ##      TARGET_TAXON,
  ## The two-letter codes for HABIT and FFG are assumed to match those used for NRSA. All names are assumed to be uppercase
  necTraits <- c("PHYLUM", "CLASS", "ORDER", "FAMILY", "SUBFAMILY", "TRIBE", "GENUS", "TARGET_TAXON", taxa_id)
  if (any(necTraits %nin% names(inTaxa))) {
    msgTraits <- which(necTraits %nin% names(inTaxa))
    return(paste("Some of the traits are missing from the taxa list. The following are \nrequired for metric calculations to run:\n", necTraits[msgTraits], "\n"))
  }

  # Rename ct and taxa_id
  names(inCts)[names(inCts) == ct] <- "TOTAL"
  names(inCts)[names(inCts) == taxa_id] <- "TAXA_ID"
  names(inTaxa)[names(inTaxa) == taxa_id] <- "TAXA_ID"

  inTaxa.1 <- subset(inTaxa, select = c("TAXA_ID", "TARGET_TAXON", "PHYLUM", "CLASS", "ORDER", "FAMILY", "SUBFAMILY", "TRIBE", "GENUS"))

  # Must first create input dataset using WSA taxonomy and traits
  inCts.1 <- merge(inCts, inTaxa.1, by = "TAXA_ID")
  inCts.1$TOTAL <- as.numeric(inCts.1$TOTAL)
  inCts.1 <- subset(inCts.1, TOTAL > 0)
  # Roll up taxa to unambiguous level for NLA
  inCts.1$TARGET_TAXON[inCts.1$TARGET_TAXON %in% c("BEZZIA/PALPOMYIA", "PROBEZZIA", "SERROMYIA", "SPHAEROMIAS", "STILOBEZZIA")] <- "CERATOPOGONINAE"
  inCts.1$TARGET_TAXON[inCts.1$TARGET_TAXON %in% c("COENAGRION/ENALLAGMA")] <- "COENAGRIONIDAE"
  inCts.1$TARGET_TAXON[inCts.1$TARGET_TAXON %in% c("CALOPARYPHUS/EUPARYPHUS")] <- "STRATIOMYIDAE"
  inCts.1$TARGET_TAXON[inCts.1$TARGET_TAXON %in% c("CHELIFERA/METACHELA")] <- "EMPIDIDAE"
  inCts.1$TARGET_TAXON[inCts.1$TARGET_TAXON %in% c("LIBELLULIDAE/CORDULIIDAE")] <- "ODONATA"
  inCts.1$TARGET_TAXON[inCts.1$TARGET_TAXON %in% c("PERICOMA/TELMATOSCOPUS")] <- "PSYCHODIDAE"
  inCts.1$TARGET_TAXON[inCts.1$TARGET_TAXON %in% c("THIENEMANNIMYIA")] <- "THIENEMANNIMYIA GENUS GR."

  # After renaming taxon, remerge with taxalist by TARGET_TAXON and use new TAXA_ID as correct
  inCts.2 <- merge(inCts.1, subset(inTaxa.1, select = c("TAXA_ID", "TARGET_TAXON")), by = "TARGET_TAXON", all.x = TRUE)
  names(inCts.2)[names(inCts.2) == "TAXA_ID.y"] <- "TAXA_ID"

  # Now sum by TAXA_ID and sample ID in case rolled up taxon names already occur in samples
  totals <- aggregate(
    x = list(TOTAL = inCts.2$TOTAL), by = inCts.2[c(sampID, "TAXA_ID")],
    FUN = sum
  )

  inCts.3 <- merge(totals, inTaxa.1, by = "TAXA_ID")
  inCts.4 <- assignDistinct(inCts.3, c(sampID),
    taxlevels = c("PHYLUM", "CLASS", "ORDER", "FAMILY", "SUBFAMILY", "TRIBE", "GENUS"),
    final.name = "TARGET_TAXON",
    special.taxa = c("THIENEMANNIMYIA GENUS GR.")
  )
  inCts.4$IS_DISTINCT <- with(inCts.4, ifelse(is.na(IS_DISTINCT), 0, IS_DISTINCT))

  outCts <- subset(inCts.4, select = c(sampID, "TAXA_ID", "TOTAL", "IS_DISTINCT"))

  return(outCts)
}
