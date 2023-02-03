#' @export
#' @title Calculate richness-related metrics
#' @description This function calculates all of the variations of
#' taxa richness that include species-level, genus-level, and
#' family-level for total counts and 300-count subsamples, for
#' both the full dataset and native taxa only. Only distinct taxa
#' are included in metrics.
#' @param indata Input data frame containing variables as identified
#' in the arguments for \emph{sampID} and \emph{inputNames}.
#' Large and rare taxa should be included in the input data.
#' @param sampID sampID A character vector containing the names of all
#' variables in indata that specify a unique sample. If not specified,
#' the default is \emph{UID}.
#' @param distVars A string with the name of the distinctness
#' variables to use for richness, which can have valid values of
#' 0 (non-distinct) or 1 (distinct).
#' @param nonnative A string with the name of the numeric variable
#' indicating a non-native taxon. A value of 1 indicates that a
#' taxon is non-native, and 0 native. All other values will be ignored.
#' @param inTaxa A string with the name of the data frame containing
#' the taxonomic data.
#' @param taxa_id A string with the name of the taxon ID variable
#' in \emph{indata} that matches that in \emph{inTaxa}. The default
#' value is \emph{TAXA_ID}.
#' @param genus A string with the name of the variable containing
#' genus name in the \emph{intaxa} data frame. Blank or missing
#' genus names will be dropped.
#' @param family A string with the name of the variable containing
#' family name in the \emph{intaxa} data frame. Blank or missing
#' family names will be dropped.
#' @param prefix A vector of the same length as distVars, containing
#' metric name modifiers for each input variable.
#'
#' @return A data frame containing the richness metrics that include
#' number of taxa for full data (native only and all taxa) and
#' 300 count subsamples (native only and all taxa), with variations
#' based on the lowest taxon possible, genus level, and family level.
#'
#' @author Karen Blocksom \email{Blocksom.Karen@epa.gov}
#'
calcZoopRichnessMetrics <- function(indata, sampID, distVars,
                           nonnative, inTaxa, taxa_id='TAXA_ID',
                           genus, family, prefix = ''){

  necVars <- c(sampID, distVars, taxa_id, nonnative)
  if(any(necVars %nin% names(indata))){
    msgTraits <- which(necVars %nin% names(indata))
    print(paste("Missing variables in input data frame:",
                paste(necVars[msgTraits], collapse=',')))
    return(NULL)
  }else{
    indata[, c(distVars, nonnative)] <- lapply(indata[, c(distVars, nonnative)], as.numeric)
  }

  necTaxVars <- c(taxa_id, genus, family)
  if(any(necTaxVars %nin% names(inTaxa))){
    msgTraits <- which(necTaxVars %nin% names(inTaxa))
    print(paste("Missing variables in input taxalist:",
                paste(necTaxVars[msgTraits], collapse=',')))
    return(NULL)
  }

  if(length(distVars) != length(prefix)){
    print("The number of values in distVars argument must be the same as in the prefix argument.")
    return(NULL)
  }

  inTaxa.1 <- subset(inTaxa, select = c(taxa_id, genus, family))

  indata.taxa <- merge(indata, inTaxa.1, by = taxa_id)

  if(length(sampID)==1){
    outdata <- data.frame(col1 = unique(indata[, sampID]))
    colnames(outdata) <- sampID
  }else{
    outdata <- as.data.frame(unique(indata[, sampID]))
  }

  for(i in 1:length(distVars)){
    ntax.full <- aggregate(x = list(TOTL_NTAX = indata[, distVars[i]]),
                      by = indata[sampID],
                      FUN = function(x){sum(x, na.rm=TRUE)})

    indata.nat <- subset(indata, eval(as.name(nonnative)) == 0)

    ntax.full.nat <- aggregate(x = list(TOTL_NAT_NTAX = indata.nat[, distVars[i]]),
                           by = indata.nat[sampID],
                           FUN = function(x){sum(x, na.rm=TRUE)})

    if(prefix[i] != ''){
      names(ntax.full) <- gsub('TOTL_', paste0('TOTL', prefix[i], '_'), names(ntax.full))
      names(ntax.full.nat) <- gsub('TOTL_', paste0('TOTL', prefix[i], '_'), names(ntax.full.nat))
    }
    # Subset to only valid GENUS values and only unique values for genus for each sample
    indata.gen <- unique(subset(indata.taxa, eval(as.name(genus))!='' & !is.na(eval(as.name(genus)))
                                & eval(as.name(distVars[i]))==1,
                         select = c(sampID, distVars[i], genus)))

    indata.gen.nat <- unique(subset(indata.taxa, eval(as.name(genus))!='' & !is.na(eval(as.name(genus))) &
                                      eval(as.name(nonnative)) == 0 & eval(as.name(distVars[i]))==1,
                             select = c(sampID, distVars[i], genus)))

    ntax.gen <- aggregate(x = list(GEN_NTAX = indata.gen[, distVars[i]]),
                          by = indata.gen[sampID],
                          FUN = function(x){sum(x, na.rm=TRUE)})

    ntax.gen.nat <- aggregate(x = list(GEN_NAT_NTAX = indata.gen.nat[, distVars[i]]),
                              by = indata.gen.nat[sampID],
                              FUN = function(x){sum(x, na.rm=TRUE)})

    if(prefix[i] != ''){
      names(ntax.gen) <- gsub('GEN_', paste0('GEN', prefix[i], '_'), names(ntax.gen))
      names(ntax.gen.nat) <- gsub('GEN_', paste0('GEN', prefix[i], '_'), names(ntax.gen.nat))
    }

    # Subset to only valid FAMILY values and only unique values for family for each sample
    indata.fam <- unique(subset(indata.taxa, eval(as.name(family))!='' & !is.na(eval(as.name(family)))
                                & eval(as.name(distVars[i]))==1,
                                select = c(sampID, distVars[i], family)))

    indata.fam.nat <- unique(subset(indata.taxa, eval(as.name(family))!='' & !is.na(eval(as.name(family))) &
                                      eval(as.name(nonnative)) == 0 & eval(as.name(distVars[i]))==1,
                                    select = c(sampID, distVars[i], family)))

    ntax.fam <- aggregate(x = list(FAM_NTAX = indata.fam[, distVars[i]]),
                          by = indata.fam[sampID],
                          FUN = function(x){sum(x, na.rm=TRUE)})

    ntax.fam.nat <- aggregate(x = list(FAM_NAT_NTAX = indata.fam.nat[, distVars[i]]),
                          by = indata.fam.nat[sampID],
                          FUN = function(x){sum(x, na.rm=TRUE)})

    if(prefix[i] != ''){
      names(ntax.fam) <- gsub('FAM_', paste0('FAM', prefix[i], '_'), names(ntax.fam))
      names(ntax.fam.nat) <- gsub('FAM_', paste0('FAM', prefix[i], '_'), names(ntax.fam.nat))
    }

    outdata <- merge(outdata, ntax.full, by = sampID) |>
      merge(ntax.gen, by = sampID, all.x=TRUE) |>
      merge(ntax.fam, by = sampID, all.x=TRUE) |>
      merge(ntax.gen.nat, by = sampID, all.x=TRUE) |>
      merge(ntax.fam.nat, by = sampID, all.x=TRUE)

  }

  var.long <- names(outdata)[!(names(outdata) %in% c(sampID))]

  outdata.long <- reshape(outdata, idvar = sampID, direction = 'long',
                                varying = var.long, timevar = 'PARAMETER',
                                v.names = 'RESULT', times = var.long)

  outdata.long$RESULT <- ifelse(is.na(outdata.long$RESULT), 0, outdata.long$RESULT)

  return(outdata.long)

}
