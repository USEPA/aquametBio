#' @export
#'
#' @title Calculate benthic MMI metrics used in NRSA MMIs
#' @description This function calculates only the benthic
#' metrics in the corresponding MMI used in the National
#' Rivers and Streams Assessment (NRSA), based on the
#' aggregated Omernik ecoregions included in the input
#' data frame.
#' @param inCts A data frame containing, at minimum, the variables
#' specified in the arguments for sampID, dist, ct, taxa_id, and
#' ecoreg. It is assumed that the data have been aggregated
#' to the taxonomic levels used in WSA/NRSA already. This can
#' be done using the function \emph{calcNRSA_BentMMImets()}.
#' @param inTaxa a data frame containing taxonomic information,
#' including variables for PHYLUM, CLASS, ORDER, FAMILY, SUBFAMILY,
#' and TRIBE, as well as autecology traits with names that match those
#' in the arguments ffg, habit, and ptv. In addition, there
#' should be a variable with the name in argument taxa_id that matches
#' with all of those in the inCts data frame. The default is the
#' bentTaxa_nrsa data frame included in the package.
#' @param sampID A character vector containing the names of all
#' variables in inCts that specify a unique sample. If not specified,
#' the default is \emph{UID}
#' @param dist A string with the name of the distinctness variable,
#' which is assumed to have only values of 0 or 1. If not specified,
#' the default is \emph{IS_DISTINCT}.
#' @param ct A string with the name of the count variable. If not
#' specified, the default is \emph{TOTAL}.
#' @param ecoreg A string with the name of the ecoregion variable.
#' Valid values that correspond to regions used in NRSA are
#' CPL, NAP, NPL, SAP, sPL, TPL, UMW, WMT, and XER.
#' @param taxa_id A string with the name of the taxon ID variable
#' in \emph{inCts} that matches that in \emph{inTaxa}. The default
#' value is \emph{TAXA_ID}.
#' @param ffg A string with the name of the functional feeding group
#' variable in inTaxa. The default value is \emph{FFG}. Values used
#' in calculations include CF, CG, PR, SH, Sc, representing
#' collector-filterer, collector-gatherer, predator, shredder, and
#' scraper, respectively. Each taxon may have more than
#' one FFG value.
#' @param habit A string with the name of the habit variable in inTaxa.
#' The default value is \emph{HABIT}. Values for habit that are used in
#' calculations include BU, CB, CN, SP, SW, representing burrower,
#' climber, clinger, sprawler, and swimmer, respectively. Each taxon
#' may have more than one value for HABIT.
#' @param ptv A string with the name of the pollution tolerance value
#' variable in inTaxa. The default is \emph{PTV}.
#' @return A data frame containing the variables in sampID and all of
#' the benthic macroinvertebrate metrics used in the MMI as additional variables.
#' The metrics generated, by aggregated ecoregion, are:
#'
#'  CPL: NOINPIND, HPRIME, SHRDNTAX, CLNGPTAX, EPT_NTAX, TOLRPTAX
#'
#'  NAP: EPT_PTAX, DOM5PIND, SCRPNTAX, CLNGPTAX, EPT_NTAX, NTOLPTAX
#'
#'  NPL: EPT_PTAX, HPRIME, SCRPNTAX, BURRPTAX, EPHENTAX, NTOLNTAX
#'
#'  SAP: EPHEPTAX, HPRIME, SCRPNTAX, BURRPTAX, EPT_NTAX, TOLRPTAX
#'
#'  SPL: EPT_PIND, HPRIME, SCRPNTAX, BURRPTAX, EPT_NTAX, INTLNTAX
#'
#'  TPL: EPT_PIND, HPRIME, SCRPNTAX, CLNGNTAX, EPHENTAX, STOLPTAX
#'
#'  UMW: CHIRPTAX, HPRIME, SHRDNTAX, BURRPTAX, EPT_NTAX, STOLPTAX
#'
#'  WMT: EPT_PTAX, DOM5PIND, SCRPNTAX, CLNGPTAX, EPT_NTAX, TOLRPTAX
#'
#'  XER: NOINPIND, DOM5PIND, SCRPNTAX, CLNGPTAX, EPT_NTAX, TOLRPTAX
#'
#' Metric descriptions are included in \emph{NRSA_Fish_Metric_Descriptions.pdf},
#' included in this package.
#' @author Karen Blocksom \email{Blocksom.Karen@epa.gov}
#' @keywords survey
calcNRSA_BentMMImets <- function(
    inCts, inTaxa = bentTaxa_nrsa, sampID = "UID", ecoreg = NULL,
    dist = "IS_DISTINCT", ct = "TOTAL", taxa_id = "TAXA_ID",
    ffg = "FFG", habit = "HABIT", ptv = "PTV") {
  inTaxa <- as.data.frame(inTaxa)
  inCts <- as.data.frame(inCts)
  # Make sure all taxa match to taxalist and send error if not
  checkTaxa <- inCts[!(inCts$TAXA_ID %in% inTaxa$TAXA_ID), ]
  if (nrow(checkTaxa) > 0) {
    return(print("Taxa in counts that do not have matches in taxalist! Cannot continue."))
  }

  if ("NON_TARGET" %in% names(inTaxa)) {
    inTaxa <- subset(inTaxa, is.na(NON_TARGET) | NON_TARGET == "" | NON_TARGET == "N")
  }
  # inTaxa <- subset(inTaxa, is.na(NON_TARGET) | NON_TARGET == "")
  inTaxa[, c(ptv, taxa_id)] <- lapply(inTaxa[, c(ptv, taxa_id)], as.numeric)

  ctVars <- c(sampID, dist, ct, taxa_id, ecoreg)
  if (any(ctVars %nin% names(inCts))) {
    msgTraits <- which(ctVars %nin% names(inCts))
    print(paste(
      "Missing variables in input data frame:",
      paste(ctVars[msgTraits], collapse = ",")
    ))
    return(NULL)
  }

  ecoCk <- unique(inCts[, ecoreg])
  ecos <- c("CPL", "NAP", "NPL", "SAP", "SPL", "TPL", "UMW", "WMT", "XER")
  if (any(ecoCk %nin% ecos)) {
    msgEco <- which(ecoCk %nin% ecos)
    print(paste(
      "These ecoregions are not valid: ",
      paste(ecoCk[msgEco], collapse = ",")
    ))
    return(NULL)
  }

  metnames <- data.frame(
    ECO = c(
      rep("CPL", 6), rep("NAP", 6), rep("NPL", 6), rep("SAP", 6), rep("SPL", 6), rep("TPL", 6),
      rep("UMW", 6), rep("WMT", 6), rep("XER", 6)
    ),
    METRIC = c(
      "NOINPIND", "HPRIME", "SHRDNTAX", "CLNGPTAX", "EPT_NTAX", "TOLRPTAX",
      "EPT_PTAX", "DOM5PIND", "SCRPNTAX", "CLNGPTAX", "EPT_NTAX", "NTOLPTAX",
      "EPT_PTAX", "HPRIME", "SCRPNTAX", "BURRPTAX", "EPHENTAX", "NTOLNTAX",
      "EPHEPTAX", "HPRIME", "SCRPNTAX", "BURRPTAX", "EPT_NTAX", "TOLRPTAX",
      "EPT_PIND", "HPRIME", "SCRPNTAX", "BURRPTAX", "EPT_NTAX", "INTLNTAX",
      "EPT_PIND", "HPRIME", "SCRPNTAX", "CLNGNTAX", "EPHENTAX", "STOLPTAX",
      "CHIRPTAX", "HPRIME", "SHRDNTAX", "BURRPTAX", "EPT_NTAX", "STOLPTAX",
      "EPT_PTAX", "DOM5PIND", "SCRPNTAX", "CLNGPTAX", "EPT_NTAX", "TOLRPTAX",
      "NOINPIND", "DOM5PIND", "SCRPNTAX", "CLNGPTAX", "EPT_NTAX", "TOLRPTAX"
    ),
    METTYPE = rep(c("TAX", "DOM", "FFG", "HAB", "TAX", "TOL"), 9),
    stringsAsFactors = FALSE
  )

  # Calculate all metrics associated with any ecoregion, then only keep those that
  # match the appropriate ecoregion for each site
  inCts <- subset(inCts, select = c(sampID, ct, dist, taxa_id, ecoreg))
  # Rename ct and dist to TOTAL and IS_DISTINCT
  names(inCts)[names(inCts) == ct] <- "TOTAL"
  names(inCts)[names(inCts) == dist] <- "IS_DISTINCT"
  names(inCts)[names(inCts) == taxa_id] <- "TAXA_ID"

  for (i in 1:length(sampID)) {
    if (i == 1) {
      inCts$SAMPID <- inCts[, sampID[i]]
    } else {
      inCts$SAMPID <- paste(inCts$SAMPID, inCts[, sampID[i]], sep = ".")
    }
  }

  samples <- unique(subset(inCts, select = c(sampID, "SAMPID", ecoreg)))

  # Taxonomy and traits checks
  necTraits <- c(
    "PHYLUM", "CLASS", "ORDER", "FAMILY", "TRIBE", "SUBFAMILY", "GENUS",
    ffg, habit, ptv
  )
  if (any(necTraits %nin% names(inTaxa))) {
    msgTraits <- which(necTraits %nin% names(inTaxa))
    return(paste("Some of the traits are missing from the taxa list. The following are \nrequired for metric calculations to run:\n", necTraits[msgTraits], "\n"))
  }

  inTaxa <- subset(inTaxa, select = names(inTaxa) %in% c(
    "TAXA_ID", "PHYLUM", "CLASS", "ORDER", "FAMILY",
    "TRIBE", "SUBFAMILY", "GENUS", ffg, habit, ptv
  ))
  names(inTaxa)[names(inTaxa) == habit] <- "HABIT"
  names(inTaxa)[names(inTaxa) == ffg] <- "FFG"
  names(inTaxa)[names(inTaxa) == ptv] <- "PTV"

  inCts.1 <- inCts[inCts$TAXA_ID %in% inTaxa$TAXA_ID, c("SAMPID", "TAXA_ID", "TOTAL", "IS_DISTINCT")]
  inCts.1 <- inCts.1[!is.na(inCts.1$TOTAL) & inCts.1$TOTAL > 0, ]

  inTaxa.1 <- inTaxa
  inTaxa.1$EPT_ <- with(inTaxa.1, ifelse(ORDER %in% c("PLECOPTERA", "EPHEMEROPTERA", "TRICHOPTERA"), 1, NA))
  inTaxa.1$EPHE <- with(inTaxa.1, ifelse(ORDER %in% c("EPHEMEROPTERA"), 1, NA))
  inTaxa.1$CHIR <- with(inTaxa.1, ifelse(FAMILY %in% c("CHIRONOMIDAE"), 1, NA))
  inTaxa.1$NOIN <- with(inTaxa.1, ifelse(CLASS %nin% c("INSECTA"), 1, NA))
  inTaxa.1$SHRD <- with(inTaxa.1, ifelse(grepl("SH", FFG), 1, NA))
  inTaxa.1$SCRP <- with(inTaxa.1, ifelse(grepl("SC", FFG), 1, NA))
  inTaxa.1$BURR <- with(inTaxa.1, ifelse(grepl("BU", HABIT), 1, NA))
  inTaxa.1$CLNG <- with(inTaxa.1, ifelse(grepl("CN", HABIT), 1, NA))
  inTaxa.1$TOLR <- with(inTaxa.1, ifelse(PTV >= 7, 1, NA))
  inTaxa.1$INTL <- with(inTaxa.1, ifelse(PTV <= 3, 1, NA))
  inTaxa.1$NTOL <- with(inTaxa.1, ifelse(PTV < 6, 1, NA))
  inTaxa.1$STOL <- with(inTaxa.1, ifelse(PTV >= 8, 1, NA))

  # Drop non-target taxa if included in taxalist
  if (length(grep("NON_TARGET", names(inTaxa.1))) > 0) {
    inTaxa.1 <- subset(inTaxa.1, is.na(NON_TARGET) | NON_TARGET == "")
  }

  params <- c("EPT_", "EPHE", "CHIR", "NOIN", "SCRP", "SHRD", "BURR", "CLNG", "TOLR", "INTL", "NTOL", "STOL")

  taxalong <- reshape(inTaxa.1[, c("TAXA_ID", params)],
    idvar = "TAXA_ID", direction = "long",
    varying = params, timevar = "TRAIT", v.names = "value",
    times = params
  )

  taxalong <- taxalong[!is.na(taxalong$value) & taxalong$value != "", ]

  taxalong$TRAIT <- as.character(taxalong$TRAIT)

  totals <- aggregate(
    x = list(TOTLNIND = inCts.1$TOTAL, TOTLNTAX = inCts.1$IS_DISTINCT), by = inCts.1[c("SAMPID")],
    FUN = sum
  )

  inCts.1 <- merge(inCts.1, totals, by = "SAMPID")
  inCts.1$CALCPIND <- with(inCts.1, TOTAL / TOTLNIND)
  inCts.1$CALCPTAX <- with(inCts.1, IS_DISTINCT / TOTLNTAX)

  # Merge the count data with the taxalist containing only the traits of
  # interest
  traitDF <- merge(inCts.1, taxalong, by = "TAXA_ID")

  # Calculate no. individuals, % individuals, no. taxa, and % taxa for each
  # trait in taxalist
  outMet.1 <- aggregate(
    x = list(NIND = traitDF$TOTAL, NTAX = traitDF$IS_DISTINCT),
    by = traitDF[c("SAMPID", "TRAIT", "TOTLNTAX")],
    FUN = sum
  )
  outMet.2 <- aggregate(
    x = list(PIND = traitDF$CALCPIND, PTAX = traitDF$CALCPTAX),
    by = traitDF[c("SAMPID", "TRAIT", "TOTLNTAX")],
    FUN = function(x) {
      round(sum(x) * 100, 2)
    }
  )

  outMet <- merge(outMet.1, outMet.2, by = c("SAMPID", "TRAIT", "TOTLNTAX"))

  outLong <- reshape(outMet,
    idvar = c("SAMPID", "TOTLNTAX", "TRAIT"), direction = "long",
    varying = names(outMet)[!names(outMet) %in% c("SAMPID", "TOTLNTAX", "TRAIT")],
    timevar = "variable", v.names = "value",
    times = names(outMet)[!names(outMet) %in% c("SAMPID", "TOTLNTAX", "TRAIT")]
  )

  outLong$variable <- paste(outLong$TRAIT, outLong$variable, sep = "")
  outLong$TRAIT <- NULL

  outWide <- reshape(outLong,
    idvar = c("SAMPID", "TOTLNTAX"), direction = "wide",
    timevar = "variable", v.names = "value"
  )

  names(outWide) <- gsub("value\\.", "", names(outWide))

  outWide <- merge(outWide, samples, by = "SAMPID")

  shanMet <- ShanDiversity(inCts.1)
  domMet <- Dominance(inCts.1, 5)
  domMet$DOM5PIND <- with(domMet, ifelse(is.na(DOM5PIND), 100, DOM5PIND))

  outAll <- merge(outWide, shanMet, by = "SAMPID")
  outAll <- merge(outAll, domMet, by = "SAMPID")

  outLong.1 <- reshape(outAll,
    idvar = c(sampID, "SAMPID", ecoreg), direction = "long",
    varying = names(outAll)[!names(outAll) %in% c(sampID, "SAMPID", ecoreg)],
    timevar = "variable", v.names = "value",
    times = names(outAll)[!names(outAll) %in% c(sampID, "SAMPID", ecoreg)]
  )

  outLong.2 <- merge(outLong.1, metnames, by.x = c(ecoreg, "variable"), by.y = c("ECO", "METRIC"))

  outLong.2$value[is.na(outLong.2$value)] <- 0

  ckMetnum <- as.data.frame(table(SAMPID = outLong.2$SAMPID))

  ckMetnum <- subset(ckMetnum, Freq != 6)
  if (nrow(ckMetnum) > 0) {
    print("Error in output! Wrong number of metrics per site!")
  }

  # Finally, we can recast the metrics df into wide format for output
  outLong.2$METTYPE <- NULL
  outWide.fin <- reshape(outLong.2,
    idvar = c(sampID, "SAMPID", ecoreg), direction = "wide",
    timevar = "variable", v.names = "value"
  )

  names(outWide.fin) <- gsub("value\\.", "", names(outWide.fin))

  outWide.fin <- merge(outWide.fin, totals[, c("SAMPID", "TOTLNIND")], by = "SAMPID")
  outWide.fin$SAMPID <- NULL

  return(outWide.fin)
}
