#' @export
#' @title Calculate taxonomy metrics for benthic macroinvertebrates
#'
#' @description This function calculates all taxonomy related
#' benthic metrics.
#'
#' @param inCts A data frame containing, at minimum, the variables
#' specified in the arguments for sampID, dist, ct, and taxa_id
#' @param inTaxa a data frame containing taxonomic information,
#' including variables for PHYLUM, CLASS, ORDER, FAMILY, SUBFAMILY,
#' and TRIBE, as well as autecology traits with names that match those
#' in the arguments ffg, habit, and ptv. In addition, there
#' should be a variable with the name in argument taxa_id that matches
#' with all of those in the indf data frame
#' @param sampID A character vector containing the names of all
#' variables in indf that specify a unique sample. If not specified,
#' the default is \emph{UID}
#' @param dist A string with the name of the distinctness variable,
#' which is assumed to have only values of 0 or 1. If not specified,
#' the default is \emph{IS_DISTINCT}.
#' @param ct A string with the name of the count variable. If not
#' specified, the default is \emph{TOTAL}.
#' @param taxa_id A string with the name of the taxon ID variable
#' in \emph{inCts} that matches that in \emph{inTaxa}. The default
#' value is \emph{TAXA_ID}.
#' @return A data frame containing the variables in sampID and
#' the benthic macroinvertebrate metrics as additional variables.
#'  These metrics are named  AMPHNTAX, AMPHPIND, AMPHPTAX, CHIRNTAX,
#'  CHIRPIND, CHIRPTAX, CRUSNTAX, CRUSPIND, CRUSPTAX, DIPTNTAX,
#'  DIPTPIND, DIPTPTAX, EPHENTAX, EPHEPIND, EPHEPTAX, EPOTNTAX,
#'  EPOTPIND, EPOTPTAX, EPT_NTAX, EPT_PIND, EPT_PTAX, HEMINTAX,
#'  HEMIPIND, HEMIPTAX, MITENTAX, MITEPIND, MITEPTAX, MOLLNTAX,
#'  MOLLPIND, MOLLPTAX, NOINNTAX, NOINPIND, NOINPTAX, ODONNTAX,
#'  ODONPIND, ODONPTAX, OLLENTAX, OLLEPIND, OLLEPTAX, ORTHNTAX,
#'  ORTHPIND, ORTHPTAX, PLECNTAX, PLECPIND, PLECPTAX, TANYNTAX,
#'  TANYPIND, TANYPTAX, TRICNTAX, TRICPIND, TRICPTAX, TUBINAIDNTAX,
#'  TUBINAIDPIND, TUBINAIDPTAX, and ORTHCHIRPIND.
#'
#' Metric descriptions are included in \emph{NRSA_Invertebrate_Metric_Descriptions.pdf},
#' included in this package.
#' @author Karen Blocksom \email{Blocksom.Karen@epa.gov}


calcBentTaxMets <- function(inCts, inTaxa, sampID="UID", dist="IS_DISTINCT",
                             ct="TOTAL",taxa_id='TAXA_ID'){

  ctVars <- c(sampID,dist,ct,taxa_id)
  if(any(ctVars %nin% names(inCts))){
    msgTraits <- which(ctVars %nin% names(inCts))
    print(paste("Missing variables in input data frame:",paste(names(inCts)[msgTraits],collapse=',')))
    return(NULL)
  }


  inCts <- subset(inCts,select=c(sampID,ct,dist,taxa_id))
  # Rename ct and dist to TOTAL and IS_DISTINCT
  names(inCts)[names(inCts)==ct] <- 'TOTAL'
  names(inCts)[names(inCts)==dist] <- 'IS_DISTINCT'
  names(inCts)[names(inCts)==taxa_id] <- 'TAXA_ID'

  for(i in 1:length(sampID)){
    if(i==1) inCts$SAMPID <- inCts[,sampID[i]]
    else inCts$SAMPID <- paste(inCts$SAMPID,inCts[,sampID[i]],sep='.')
  }

  # Make sure all taxa match to taxalist and send error if not
  checkTaxa <- inCts[!(inCts$TAXA_ID %in% inTaxa$TAXA_ID),]
  if(nrow(checkTaxa)>0){
    return(print('Taxa in counts that do not have matches in taxalist! Cannot continue.'))
  }

  # If necessary, load the bentTaxa data frame and assign it to inTaxa.  Though
  # NON_TARGET taxa are included in the table provided by NRSA, we need to exclude
  # them from our calculations.
  if('NON_TARGET' %in% names(inTaxa)){
    inTaxa <- subset(inTaxa, is.na(NON_TARGET) | NON_TARGET == "" |NON_TARGET=='N')
  }
  # }

  ## This code assumes that the following are columns in the taxa file: PHYLUM, CLASS, ORDER, FAMILY, SUBFAMILY, TRIBE, HABIT, FFG, PTV,
  ##    	TARGET_TAXON,
  ## The two-letter codes for HABIT and FFG are assumed to match those used for NRSA. All names are assumed to be uppercase
  necTraits <- c('PHYLUM','CLASS','ORDER','FAMILY','TRIBE','SUBFAMILY','GENUS')
  if(any(necTraits %nin% names(inTaxa))){
    msgTraits <- which(necTraits %nin% names(inTaxa))
    return(paste("Some of the traits are missing from the taxa list. The following are \nrequired for metric calculations to run:\n", necTraits[msgTraits], "\n"))
  }

  samples <- unique(subset(inCts,select=c(sampID,'SAMPID')))

  inCts.1 <- inCts[inCts$TAXA_ID %in% inTaxa$TAXA_ID, c('SAMPID','TAXA_ID','TOTAL','IS_DISTINCT')]

  inCts.1 <- inCts.1[!is.na(inCts.1$TOTAL) & inCts.1$TOTAL>0,]

  inTaxa.1 <- inTaxa
  inTaxa.1$EPT_ <- with(inTaxa.1, ifelse(ORDER %in% c('PLECOPTERA','EPHEMEROPTERA','TRICHOPTERA'),1,NA))
  inTaxa.1$EPHE <- with(inTaxa.1, ifelse(ORDER %in% c('EPHEMEROPTERA'),1,NA))
  inTaxa.1$PLEC <- with(inTaxa.1, ifelse(ORDER %in% c('PLECOPTERA'),1,NA))
  inTaxa.1$TRIC <- with(inTaxa.1, ifelse(ORDER %in% c('TRICHOPTERA'),1,NA))
  inTaxa.1$CHIR <- with(inTaxa.1, ifelse(FAMILY %in% c('CHIRONOMIDAE'),1,NA))
  inTaxa.1$CRUS <- with(inTaxa.1, ifelse(CLASS %in% c('MALACOSTRACA','MAXILLOPODA','BRANCHIOPODA'
                            ,'CEPHALOCARIDA','OSTRACODA','REMIPEDIA'),1,NA))
  inTaxa.1$NOIN <- with(inTaxa.1, ifelse(CLASS %nin% c('INSECTA'),1,NA))
  inTaxa.1$DIPT <- with(inTaxa.1, ifelse(ORDER %in% c('DIPTERA'),1,NA))
  inTaxa.1$MOLL <- with(inTaxa.1, ifelse(PHYLUM %in% c('MOLLUSCA'),1,NA))
  inTaxa.1$AMPH <- with(inTaxa.1, ifelse(ORDER %in% c('AMPHIPODA'),1,NA))
  inTaxa.1$EPOT <- with(inTaxa.1, ifelse(ORDER %in% c('ODONATA','PLECOPTERA','EPHEMEROPTERA','TRICHOPTERA'),1,NA))
  inTaxa.1$HEMI <- with(inTaxa.1, ifelse(ORDER %in% c('HEMIPTERA'),1,NA))
  inTaxa.1$MITE <- with(inTaxa.1, ifelse(ORDER %in% c('TROMBIDIFORMES','SARCOPTIFORMES'),1,NA))
  inTaxa.1$ODON <- with(inTaxa.1, ifelse(ORDER %in% c('ODONATA'),1,NA))
  inTaxa.1$OLLE <- with(inTaxa.1, ifelse(CLASS %in% c('OLIGOCHAETA','HIRUDINEA','CLITELLATA'),1,NA))
  inTaxa.1$ORTH <- with(inTaxa.1, ifelse(SUBFAMILY %in% c('ORTHOCLADIINAE'),1,NA))
  inTaxa.1$TANY <- with(inTaxa.1, ifelse(TRIBE %in% c('TANYTARSINI'),1,NA))
  inTaxa.1$TUBINAID <- with(inTaxa.1, ifelse(FAMILY %in% c('TUBIFICIDAE','NAIDIDAE'),1,NA))

  # Drop non-target taxa if included in taxalist
  if(length(grep('NON_TARGET',names(inTaxa.1)))>0) {
    inTaxa.1 <- subset(inTaxa.1,is.na(NON_TARGET)|NON_TARGET=='')
  }

  params<-c('EPT_','EPHE','PLEC','TRIC','CHIR','CRUS','DIPT','MOLL','NOIN'
            ,'TUBINAID','AMPH','EPOT','HEMI','MITE','ODON','OLLE','ORTH','TANY')


  taxalong <- reshape(inTaxa.1[,c('TAXA_ID',params)], idvar = 'TAXA_ID', direction = 'long',
                      varying = params, timevar = 'TRAIT', v.names = 'value',
                      times = params)

  taxalong <- taxalong[!is.na(taxalong$value),]

  taxalong$TRAIT <- as.character(taxalong$TRAIT)

  totals <- aggregate(x = list(TOTLNIND = inCts.1$TOTAL, TOTLNTAX = inCts.1$IS_DISTINCT), by = inCts.1[c('SAMPID')],
                       FUN = sum)

  inCts.1 <- merge(inCts.1, totals, by = 'SAMPID')
  inCts.1$CALCPIND <- with(inCts.1, TOTAL/TOTLNIND)
  inCts.1$CALCPTAX <- with(inCts.1, IS_DISTINCT/TOTLNTAX)

  # Merge the count data with the taxalist containing only the traits of
  # interest
  traitDF <- merge(inCts.1, taxalong, by='TAXA_ID')

  # Calculate no. individuals, % individuals, no. taxa, and % taxa for each
  # trait in taxalist
  outMet.1 <- aggregate(x = list(NIND = traitDF$TOTAL, NTAX = traitDF$IS_DISTINCT),
                                    by = traitDF[c('SAMPID','TRAIT','TOTLNTAX')],
                      FUN = sum)
  outMet.2 <- aggregate(x = list(PIND = traitDF$CALCPIND, PTAX = traitDF$CALCPTAX),
                        by = traitDF[c('SAMPID','TRAIT','TOTLNTAX')],
                        FUN = function(x){round(sum(x)*100, 2)})

  outMet <- merge(outMet.1, outMet.2, by = c('SAMPID','TRAIT','TOTLNTAX'))

  # Melt df to create metric names, then recast into wide format with metric
  # names
  print("Done calculating taxonomy metrics.")

  outLong <- reshape(outMet, idvar = c('SAMPID','TOTLNTAX','TRAIT'), direction = 'long',
                     varying = names(outMet)[!names(outMet) %in% c('SAMPID','TOTLNTAX','TRAIT')],
                     timevar = 'variable', v.names = 'value',
                     times = names(outMet)[!names(outMet) %in% c('SAMPID','TOTLNTAX','TRAIT')])

  outLong$variable <- paste(outLong$TRAIT,outLong$variable,sep='')
  outLong$TRAIT <- NULL

  outWide <- reshape(outLong, idvar = c('SAMPID','TOTLNTAX'), direction = 'wide',
                     timevar = 'variable', v.names = 'value')

  names(outWide) <- gsub("value\\.", "", names(outWide))

  outWide <- merge(outWide, samples, by='SAMPID')


  # ORTHCHIRPIND are % of Chironomidae individuals in ORTHOCLADIINAE (not in
  # total indiv in sample)
  outWide$ORTHCHIRPIND <- with(outWide, round(ORTHNIND/CHIRNIND*100,2))
  outWide$ORTHCHIRPIND <- with(outWide, ifelse(is.na(ORTHCHIRPIND)|is.nan(ORTHCHIRPIND),0,ORTHCHIRPIND))

  empty_tax <- data.frame(t(rep(NA,56)),stringsAsFactors=F)
  names(empty_tax) <- c('TOTLNTAX','AMPHNTAX','AMPHPIND'
                        ,'AMPHPTAX','CHIRNTAX','CHIRPIND'
                        ,'CHIRPTAX'
                        ,'CRUSNTAX','CRUSPIND','CRUSPTAX'
                        ,'DIPTNTAX','DIPTPIND','DIPTPTAX','EPHENTAX'
                        ,'EPHEPIND','EPHEPTAX','EPOTNTAX','EPOTPIND','EPOTPTAX','EPT_NTAX','EPT_PIND'
                        ,'EPT_PTAX','HEMINTAX','HEMIPIND','HEMIPTAX'
                        ,'MITENTAX','MITEPIND','MITEPTAX'
                        ,'MOLLNTAX','MOLLPIND','MOLLPTAX','NOINNTAX','NOINPIND','NOINPTAX'
                        ,'ODONNTAX','ODONPIND','ODONPTAX','OLLENTAX','OLLEPIND'
                        ,'OLLEPTAX','ORTHNTAX','ORTHPIND','ORTHPTAX','PLECNTAX','PLECPIND','PLECPTAX'
                        ,'TANYNTAX','TANYPIND'
                        ,'TANYPTAX','TRICNTAX','TRICPIND','TRICPTAX'
                        ,'TUBINAIDNTAX','TUBINAIDPIND','TUBINAIDPTAX','ORTHCHIRPIND')

  outWide.all <- merge(outWide, empty_tax, all=TRUE)
  outWide.all <- subset(outWide.all, !is.na(SAMPID)) # This step also drops variables with all missing values (HEMI metrics)

  # If we re-melt df now, we have missing values where the metric should be a
  # zero, so we can set NAs to 0 now
  outWide.all[is.na(outWide.all)] <- 0

  outWide.fin <- outWide.all
  outWide.fin$SAMPID <- NULL

  outWide.fin <- outWide.fin[, c(sampID,'AMPHNTAX','AMPHPIND'
                                 ,'AMPHPTAX','CHIRNTAX','CHIRPIND'
                                 ,'CHIRPTAX'
                                 ,'CRUSNTAX','CRUSPIND','CRUSPTAX'
                                 ,'DIPTNTAX','DIPTPIND','DIPTPTAX','EPHENTAX'
                                 ,'EPHEPIND','EPHEPTAX','EPOTNTAX','EPOTPIND','EPOTPTAX','EPT_NTAX','EPT_PIND'
                                 ,'EPT_PTAX','HEMINTAX','HEMIPIND','HEMIPTAX'
                                 ,'MITENTAX','MITEPIND','MITEPTAX'
                                 ,'MOLLNTAX','MOLLPIND','MOLLPTAX','NOINNTAX','NOINPIND','NOINPTAX'
                                 ,'ODONNTAX','ODONPIND','ODONPTAX','OLLENTAX','OLLEPIND'
                                 ,'OLLEPTAX','ORTHNTAX','ORTHPIND','ORTHPTAX','PLECNTAX','PLECPIND','PLECPTAX'
                                 ,'TANYNTAX','TANYPIND'
                                 ,'TANYPTAX','TRICNTAX','TRICPIND','TRICPTAX'
                                 ,'TUBINAIDNTAX','TUBINAIDPIND','TUBINAIDPTAX','ORTHCHIRPIND')]

  return(outWide.fin)

}


#' @export
#' @title Calculate functional feeding group benthic metrics
#'
#' @description This function calculates
#' all functional feeding group benthic metrics using trait
#' information supplied.
#'
#' @param inCts A data frame containing, at minimum, the variables
#' specified in the arguments for sampID, dist, ct, and taxa_id
#' @param inTaxa a data frame containing taxonomic information,
#' including variables for PHYLUM, CLASS, ORDER, FAMILY, SUBFAMILY,
#' and TRIBE, as well as autecology traits with names that match those
#' in the arguments ffg, habit, and ptv. In addition, there
#' should be a variable with the name in argument taxa_id that matches
#' with all of those in the indf data frame
#' @param sampID A character vector containing the names of all
#' variables in indf that specify a unique sample. If not specified,
#' the default is \emph{UID}
#' @param dist A string with the name of the distinctness variable,
#' which is assumed to have only values of 0 or 1. If not specified,
#' the default is \emph{IS_DISTINCT}.
#' @param ct A string with the name of the count variable. If not
#' specified, the default is \emph{TOTAL}.
#' @param taxa_id A string with the name of the taxon ID variable
#' in \emph{inCts} that matches that in \emph{inTaxa}. The default
#' value is \emph{TAXA_ID}.
#' @param ffg A string with the name of the functional feeding group
#' variable in inTaxa. The default value is \emph{FFG}. Values used
#' in calculations include CF, CG, PR, SH, Sc, representing
#' collector-filterer, collector-gatherer, predator, shredder, and
#' scraper, respectively. Each taxon may have more than
#' one ffg value.
#' @return A data frame containing the variables in sampID and
#' the benthic macroinvertebrate metrics as additional variables.
#' Metrics included are named COFINTAX, COFIPIND, COFIPTAX,
#' COFITRICNTAX, COFITRICPIND, COFITRICPTAX, COGANTAX, COGAPIND,
#' COGAPTAX, PREDNTAX, PREDPIND, PREDPTAX, SCRPNTAX, SCRPPIND,
#' SCRPPTAX, SHRDNTAX, SHRDPIND, and SHRDPTAX.
#'
#' Metric descriptions are included in \emph{NRSA_Invertebrate_Metric_Descriptions.pdf},
#' included in this package.
#' @author Karen Blocksom \email{Blocksom.Karen@epa.gov}

calcBentFFGmets <- function(inCts, inTaxa, sampID="UID", dist="IS_DISTINCT",
                        ct="TOTAL",taxa_id='TAXA_ID',ffg='FFG'){
  ctVars <- c(sampID,dist,ct,taxa_id)
  if(any(ctVars %nin% names(inCts))){
    msgTraits <- which(ctVars %nin% names(inCts))
    print(paste("Missing variables in input data frame:",paste(names(inCts)[msgTraits],collapse=',')))
    return(NULL)
  }

  inCts <- subset(inCts,select=c(sampID,ct,dist,taxa_id))
  # Rename ct and dist to TOTAL and IS_DISTINCT
  names(inCts)[names(inCts)==ct] <- 'TOTAL'
  names(inCts)[names(inCts)==dist] <- 'IS_DISTINCT'
  names(inCts)[names(inCts)==taxa_id] <- 'TAXA_ID'

  for(i in 1:length(sampID)){
    if(i==1) inCts$SAMPID <- inCts[,sampID[i]]
    else inCts$SAMPID <- paste(inCts$SAMPID,inCts[,sampID[i]],sep='.')
  }

  # Make sure all taxa match to taxalist and send error if not
  checkTaxa <- inCts[!(inCts$TAXA_ID %in% inTaxa$TAXA_ID),]
  if(nrow(checkTaxa)>0){
    return(print('Taxa in counts that do not have matches in taxalist! Cannot continue.'))
  }

  # If necessary, load the bentTaxa data frame and assign it to inTaxa.  Though
  # NON_TARGET taxa are included in the table provided by NRSA, we need to exclude
  # them from our calculations.
  if('NON_TARGET' %in% names(inTaxa)){
    inTaxa <- subset(inTaxa, is.na(NON_TARGET) | NON_TARGET == "" |NON_TARGET=='N')
  }

  if(ffg %nin% names(inTaxa)){
    return(paste("Input taxa list does not contain variable",ffg,"."))
  }

  if('ORDER' %nin% names(inTaxa)){
    print("Missing variable ORDER from input taxa list. Will not calculate Collector-filterer Trichoptera metrics.")
  }

  inTaxa <- subset(inTaxa,select=names(inTaxa) %in% c('TAXA_ID',ffg,'ORDER'))
  names(inTaxa)[names(inTaxa)==ffg] <- 'FFG'

  samples <- unique(subset(inCts,select=c(sampID,'SAMPID')))

  inCts.1 <- inCts[inCts$TAXA_ID %in% inTaxa$TAXA_ID, c('SAMPID','TAXA_ID','TOTAL','IS_DISTINCT')]
  inCts.1 <- inCts.1[!is.na(inCts.1$TOTAL) & inCts.1$TOTAL>0,]

  inTaxa.1 <- inTaxa
  inTaxa.1$COFI <- with(inTaxa.1, ifelse(grepl('CF',FFG), 1, NA))
  inTaxa.1$COGA <- with(inTaxa.1, ifelse(grepl('CG',FFG), 1, NA))
  inTaxa.1$PRED <- with(inTaxa.1, ifelse(grepl('PR',FFG), 1, NA))
  inTaxa.1$SHRD <- with(inTaxa.1, ifelse(grepl('SH',FFG), 1, NA))
  inTaxa.1$SCRP <- with(inTaxa.1, ifelse(grepl('SC',FFG), 1, NA))

  if('ORDER' %in% names(inTaxa.1)){
    inTaxa.1$COFITRIC <- with(inTaxa.1, ifelse(grepl('CF',FFG) & ORDER=='TRICHOPTERA', 1, NA))
  }

  # Drop non-target taxa if included in taxalist
  if(length(grep('NON_TARGET',names(inTaxa.1)))>0) {
    inTaxa.1 <- subset(inTaxa.1,is.na(NON_TARGET)|NON_TARGET=='')
  }

  params<-c('COFI','COGA','PRED','SHRD','SCRP','COFITRIC')

  taxalong <- reshape(inTaxa.1[,c('TAXA_ID',params)], idvar = 'TAXA_ID', direction = 'long',
                      varying = params, timevar = 'TRAIT', v.names = 'value',
                      times = params)

  taxalong <- taxalong[!is.na(taxalong$value),]

  taxalong$TRAIT <- as.character(taxalong$TRAIT)

  totals <- aggregate(x = list(TOTLNIND = inCts.1$TOTAL, TOTLNTAX = inCts.1$IS_DISTINCT), by = inCts.1[c('SAMPID')],
                      FUN = sum)

  inCts.1 <- merge(inCts.1, totals, by = 'SAMPID')
  inCts.1$CALCPIND <- with(inCts.1, TOTAL/TOTLNIND)
  inCts.1$CALCPTAX <- with(inCts.1, IS_DISTINCT/TOTLNTAX)

  # Merge the count data with the taxalist containing only the traits of
  # interest
  traitDF <- merge(inCts.1, taxalong, by='TAXA_ID')

  # Calculate no. individuals, % individuals, no. taxa, and % taxa for each
  # trait in taxalist
  outMet.1 <- aggregate(x = list(NIND = traitDF$TOTAL, NTAX = traitDF$IS_DISTINCT),
                        by = traitDF[c('SAMPID','TRAIT','TOTLNTAX')],
                        FUN = sum)
  outMet.2 <- aggregate(x = list(PIND = traitDF$CALCPIND, PTAX = traitDF$CALCPTAX),
                        by = traitDF[c('SAMPID','TRAIT','TOTLNTAX')],
                        FUN = function(x){round(sum(x)*100, 2)})

  outMet <- merge(outMet.1, outMet.2, by = c('SAMPID','TRAIT','TOTLNTAX'))

  # # Calculate no. individuals, % individuals, no. taxa, and % taxa for each
  # # trait in taxalist
  # Melt df to create metric names, then recast into wide format with metric
  # names
  print("Done calculating functional feeding group metrics.")

  outLong <- reshape(outMet, idvar = c('SAMPID','TOTLNTAX','TRAIT'), direction = 'long',
                     varying = names(outMet)[!names(outMet) %in% c('SAMPID','TOTLNTAX','TRAIT')],
                     timevar = 'variable', v.names = 'value',
                     times = names(outMet)[!names(outMet) %in% c('SAMPID','TOTLNTAX','TRAIT')])

  outLong$variable <- paste(outLong$TRAIT,outLong$variable,sep='')
  outLong$TRAIT <- NULL

  outWide <- reshape(outLong, idvar = c('SAMPID','TOTLNTAX'), direction = 'wide',
                     timevar = 'variable', v.names = 'value')

  names(outWide) <- gsub("value\\.", "", names(outWide))

  outWide <- merge(outWide, samples, by='SAMPID')

  empty_tax <- data.frame(t(rep(NA,18)),stringsAsFactors=F)
  names(empty_tax) <- c('COFINTAX','COFIPIND','COFIPTAX'
                        ,'COFITRICNTAX','COFITRICPIND','COFITRICPTAX'
                        ,'COGANTAX','COGAPIND','COGAPTAX'
                        ,'PREDNTAX','PREDPIND','PREDPTAX'
                        ,'SCRPNTAX','SCRPPIND','SCRPPTAX'
                        ,'SHRDNTAX','SHRDPIND','SHRDPTAX')

  outWide.all <- merge(outWide, empty_tax, all=TRUE)
  outWide.all <- subset(outWide.all, !is.na(SAMPID))

  # # If we re-melt df now, we have missing values where the metric should be a
  # # zero, so we can set NAs to 0 now
  outWide.all[is.na(outWide.all)] <- 0
  # # Finally, we can recast the metrics df into wide format for output
  outWide.fin <- outWide.all
  outWide.fin$SAMPID <- NULL

  outWide.fin <- outWide.fin[, c(sampID,'COFINTAX','COFIPIND','COFIPTAX'
                                 ,'COFITRICNTAX','COFITRICPIND','COFITRICPTAX'
                                 ,'COGANTAX','COGAPIND','COGAPTAX'
                                 ,'PREDNTAX','PREDPIND','PREDPTAX'
                                 ,'SCRPNTAX','SCRPPIND','SCRPPTAX'
                                 ,'SHRDNTAX','SHRDPIND','SHRDPTAX')]

  return(outWide.fin)

}

#' @export
#' @title Calculate habit benthic metrics
#'
#' @description This function calculates habit
#' benthic metrics.
#'
#' @param inCts A data frame containing, at minimum, the variables
#' specified in the arguments for sampID, dist, ct, and taxa_id
#' @param inTaxa a data frame containing taxonomic information,
#' including variables for PHYLUM, CLASS, ORDER, FAMILY, SUBFAMILY,
#' and TRIBE, as well as autecology traits with names that match those
#' in the arguments ffg, habit, and ptv. In addition, there
#' should be a variable with the name in argument taxa_id that matches
#' with all of those in the indf data frame
#' @param sampID A character vector containing the names of all
#' variables in indf that specify a unique sample. If not specified,
#' the default is \emph{UID}
#' @param dist A string with the name of the distinctness variable,
#' which is assumed to have only values of 0 or 1. If not specified,
#' the default is \emph{IS_DISTINCT}.
#' @param ct A string with the name of the count variable. If not
#' specified, the default is \emph{TOTAL}.
#' @param taxa_id A string with the name of the taxon ID variable
#' in \emph{inCts} that matches that in \emph{inTaxa}. The default
#' value is \emph{TAXA_ID}.
#' @param habit A string with the name of the habit variable in inTaxa.
#' The default value is \emph{HABIT}. Values for habit that are used in
#' calculations include BU, CB, CN, SP, SW, representing burrower,
#' climber, clinger, sprawler, and swimmer, respectively. Each taxon
#' may have more than one value for habit.
#' @return A data frame containing the variables in sampID and
#' the benthic macroinvertebrate metrics as additional variables.
#'  The metric names include  BURRNTAX, BURRPIND, BURRPTAX, CLMBNTAX,
#'  CLMBPIND, CLMBPTAX, CLNGNTAX, CLNGPIND, CLNGPTAX, SPWLNTAX,
#'  SPWLPIND, SPWLPTAX, SWIMNTAX, SWIMPIND, and SWIMPTAX.
#'
#' Metric descriptions are included in \emph{NRSA_Invertebrate_Metric_Descriptions.pdf},
#' included in this package.
#' @author Karen Blocksom \email{Blocksom.Karen@epa.gov}

calcBentHabitMets <- function(inCts, inTaxa, sampID="UID", dist="IS_DISTINCT",
                          ct="TOTAL",taxa_id='TAXA_ID',habit='HABIT'){

  ctVars <- c(sampID,dist,ct,taxa_id)
  if(any(ctVars %nin% names(inCts))){
    msgTraits <- which(ctVars %nin% names(inCts))
    print(paste("Missing variables in input data frame:",paste(names(inCts)[msgTraits],collapse=',')))
    return(NULL)
  }

  inCts <- subset(inCts,select=c(sampID,ct,dist,taxa_id))
  # Rename ct and dist to TOTAL and IS_DISTINCT
  names(inCts)[names(inCts)==ct] <- 'TOTAL'
  names(inCts)[names(inCts)==dist] <- 'IS_DISTINCT'
  names(inCts)[names(inCts)==taxa_id] <- 'TAXA_ID'

  for(i in 1:length(sampID)){
    if(i==1) inCts$SAMPID <- inCts[,sampID[i]]
    else inCts$SAMPID <- paste(inCts$SAMPID,inCts[,sampID[i]],sep='.')
  }

  # Make sure all taxa match to taxalist and send error if not
  checkTaxa <- inCts[!(inCts$TAXA_ID %in% inTaxa$TAXA_ID),]
  if(nrow(checkTaxa)>0){
    return(print('Taxa in counts that do not have matches in taxalist! Cannot continue.'))
  }

  # If necessary, load the bentTaxa data frame and assign it to inTaxa.  Though
  # NON_TARGET taxa are included in the table provided by NRSA, we need to exclude
  # them from our calculations.
  if('NON_TARGET' %in% names(inTaxa)){
    inTaxa <- subset(inTaxa, is.na(NON_TARGET) | NON_TARGET == "" |NON_TARGET=='N')
  }
  # }

  if(habit %nin% names(inTaxa)){
    return(paste("Input taxa list does not contain variable",habit,"."))
  }

  inTaxa <- subset(inTaxa,select=names(inTaxa) %in% c('TAXA_ID',habit))
  names(inTaxa)[names(inTaxa)==habit] <- 'HABIT'

  samples <- unique(subset(inCts,select=c(sampID,'SAMPID')))

  inCts.1 <- inCts[inCts$TAXA_ID %in% inTaxa$TAXA_ID, c('SAMPID','TAXA_ID','TOTAL','IS_DISTINCT')]
  inCts.1 <- inCts.1[!is.na(inCts.1$TOTAL) & inCts.1$TOTAL>0,]

  inTaxa.1 <- inTaxa
  inTaxa.1$BURR <- with(inTaxa.1, ifelse(grepl('BU',HABIT), 1, NA))
  inTaxa.1$CLMB <- with(inTaxa.1, ifelse(grepl('CB',HABIT), 1, NA))
  inTaxa.1$CLNG <- with(inTaxa.1, ifelse(grepl('CN',HABIT), 1, NA))
  inTaxa.1$SPWL <- with(inTaxa.1, ifelse(grepl('SP',HABIT), 1, NA))
  inTaxa.1$SWIM <- with(inTaxa.1, ifelse(grepl('SW',HABIT), 1, NA))

  # Drop non-target taxa if included in taxalist
  if(length(grep('NON_TARGET',names(inTaxa.1)))>0) {
    inTaxa.1 <- subset(inTaxa.1,is.na(NON_TARGET)|NON_TARGET=='')
  }

  params<-c('BURR','CLMB','CLNG','SPWL','SWIM')

  taxalong <- reshape(inTaxa.1[,c('TAXA_ID',params)], idvar = 'TAXA_ID', direction = 'long',
                      varying = params, timevar = 'TRAIT', v.names = 'value',
                      times = params)

  taxalong <- taxalong[!is.na(taxalong$value),]

  taxalong$TRAIT <- as.character(taxalong$TRAIT)

  totals <- aggregate(x = list(TOTLNIND = inCts.1$TOTAL, TOTLNTAX = inCts.1$IS_DISTINCT), by = inCts.1[c('SAMPID')],
                      FUN = sum)

  inCts.1 <- merge(inCts.1, totals, by = 'SAMPID')
  inCts.1$CALCPIND <- with(inCts.1, TOTAL/TOTLNIND)
  inCts.1$CALCPTAX <- with(inCts.1, IS_DISTINCT/TOTLNTAX)

  # Merge the count data with the taxalist containing only the traits of
  # interest
  traitDF <- merge(inCts.1, taxalong, by='TAXA_ID')

  # Calculate no. individuals, % individuals, no. taxa, and % taxa for each
  # trait in taxalist
  outMet.1 <- aggregate(x = list(NIND = traitDF$TOTAL, NTAX = traitDF$IS_DISTINCT),
                        by = traitDF[c('SAMPID','TRAIT','TOTLNTAX')],
                        FUN = sum)
  outMet.2 <- aggregate(x = list(PIND = traitDF$CALCPIND, PTAX = traitDF$CALCPTAX),
                        by = traitDF[c('SAMPID','TRAIT','TOTLNTAX')],
                        FUN = function(x){round(sum(x)*100, 2)})

  outMet <- merge(outMet.1, outMet.2, by = c('SAMPID','TRAIT','TOTLNTAX'))

  # Melt df to create metric names, then recast into wide format with metric
  # names
  print("Done calculating habit metrics.")

  outLong <- reshape(outMet, idvar = c('SAMPID','TOTLNTAX','TRAIT'), direction = 'long',
                     varying = names(outMet)[!names(outMet) %in% c('SAMPID','TOTLNTAX','TRAIT')],
                     timevar = 'variable', v.names = 'value',
                     times = names(outMet)[!names(outMet) %in% c('SAMPID','TOTLNTAX','TRAIT')])

  outLong$variable <- paste(outLong$TRAIT,outLong$variable,sep='')
  outLong$TRAIT <- NULL

  outWide <- reshape(outLong, idvar = c('SAMPID','TOTLNTAX'), direction = 'wide',
                     timevar = 'variable', v.names = 'value')

  names(outWide) <- gsub("value\\.", "", names(outWide))

  outWide <- merge(outWide, samples, by='SAMPID')

  empty_tax <- data.frame(t(rep(NA,15)),stringsAsFactors=F)
  names(empty_tax) <- c('BURRNTAX','BURRPIND','BURRPTAX'
                        ,'CLMBNTAX','CLMBPIND','CLMBPTAX'
                        ,'CLNGNTAX','CLNGPIND','CLNGPTAX'
                        ,'SPWLNTAX','SPWLPIND','SPWLPTAX'
                        ,'SWIMNTAX','SWIMPIND','SWIMPTAX')

  outWide.all <- merge(outWide, empty_tax, all=TRUE)
  outWide.all <- subset(outWide.all, !is.na(SAMPID))

  # # If we re-melt df now, we have missing values where the metric should be a
  # # zero, so we can set NAs to 0 now
  outWide.all[is.na(outWide.all)] <- 0
  # # Finally, we can recast the metrics df into wide format for output
  outWide.fin <- outWide.all
  outWide.fin$SAMPID <- NULL

  outWide.fin <- outWide.fin[, c(sampID,'BURRNTAX','BURRPIND','BURRPTAX'
                                 ,'CLMBNTAX','CLMBPIND','CLMBPTAX'
                                 ,'CLNGNTAX','CLNGPIND','CLNGPTAX'
                                 ,'SPWLNTAX','SPWLPIND','SPWLPTAX'
                                 ,'SWIMNTAX','SWIMPIND','SWIMPTAX')]

  return(outWide.fin)

}

#' @export
#' @title Calculate tolerance benthic metrics
#'
#' @description This function calculates tolerance
#' benthic metrics.
#'
#' @param inCts A data frame containing, at minimum, the variables
#' specified in the arguments for sampID, dist, ct, and taxa_id
#' @param inTaxa a data frame containing taxonomic information,
#' including variables for PHYLUM, CLASS, ORDER, FAMILY, SUBFAMILY,
#' and TRIBE, as well as autecology traits with names that match those
#' in the arguments ffg, habit, and ptv. In addition, there
#' should be a variable with the name in argument taxa_id that matches
#' with all of those in the indf data frame
#' @param sampID A character vector containing the names of all
#' variables in indf that specify a unique sample. If not specified,
#' the default is \emph{UID}
#' @param dist A string with the name of the distinctness variable,
#' which is assumed to have only values of 0 or 1. If not specified,
#' the default is \emph{IS_DISTINCT}.
#' @param ct A string with the name of the count variable. If not
#' specified, the default is \emph{TOTAL}.
#' @param taxa_id A string with the name of the taxon ID variable
#' in \emph{inCts} that matches that in \emph{inTaxa}. The default
#' value is \emph{TAXA_ID}.
#' @param ptv A string with the name of the pollution tolerance value
#' variable in inTaxa. The default value is \emph{PTV}. Valid values
#' for ptv range from 0 to 10.
#' @return A data frame containing the variables in sampID and
#' the benthic macroinvertebrate metrics as additional variables.
#'  The names of metrics are  FACLNTAX, FACLPIND, FACLPTAX, INTLNTAX,
#'  INTLPIND, INTLPTAX, NTOLNTAX, NTOLPIND, NTOLPTAX, STOLNTAX,
#'  STOLPIND, STOLPTAX, TL01NTAX, TL01PIND, TL01PTAX, TL23NTAX,
#'  TL23PIND, TL23PTAX, TL45NTAX, TL45PIND, TL45PTAX, TL67NTAX,
#'  TL67PIND, TL67PTAX, TOLRNTAX, TOLRPIND, TOLRPTAX, and WTD_TV.
#'
#' Metric descriptions are included in \emph{NRSA_Invertebrate_Metric_Descriptions.pdf},
#' included in this package.
#' @author Karen Blocksom \email{Blocksom.Karen@epa.gov}

calcBentTolMets <- function(inCts, inTaxa, sampID="UID", dist="IS_DISTINCT",
                          ct="TOTAL",taxa_id='TAXA_ID',ptv='PTV'){
  ctVars <- c(sampID,dist,ct,taxa_id)
  if(any(ctVars %nin% names(inCts))){
    msgTraits <- which(ctVars %nin% names(inCts))
    print(paste("Missing variables in input data frame:",paste(names(inCts)[msgTraits],collapse=',')))
    return(NULL)
  }

  inCts <- subset(inCts,select=c(sampID,ct,dist,taxa_id))
  # Rename ct and dist to TOTAL and IS_DISTINCT
  names(inCts)[names(inCts)==ct] <- 'TOTAL'
  names(inCts)[names(inCts)==dist] <- 'IS_DISTINCT'
  names(inCts)[names(inCts)==taxa_id] <- 'TAXA_ID'

  for(i in 1:length(sampID)){
    if(i==1) inCts$SAMPID <- inCts[,sampID[i]]
    else inCts$SAMPID <- paste(inCts$SAMPID,inCts[,sampID[i]],sep='.')
  }

  # Make sure all taxa match to taxalist and send error if not
  checkTaxa <- inCts[!(inCts$TAXA_ID %in% inTaxa$TAXA_ID),]
  if(nrow(checkTaxa)>0){
    return(print('Taxa in counts that do not have matches in taxalist! Cannot continue.'))
  }

  # If necessary, load the bentTaxa data frame and assign it to inTaxa.  Though
  # NON_TARGET taxa are included in the table provided by NRSA, we need to exclude
  # them from our calculations.
  if('NON_TARGET' %in% names(inTaxa)){
    inTaxa <- subset(inTaxa, is.na(NON_TARGET) | NON_TARGET == "" |NON_TARGET=='N')
  }
  # }

  if(ptv %nin% names(inTaxa)){
    return(paste("Input taxa list does not contain variable",ptv,"."))
  }

  inTaxa <- subset(inTaxa,select=names(inTaxa) %in% c('TAXA_ID',ptv))
  names(inTaxa)[names(inTaxa)==ptv] <- 'PTV'
  inTaxa$PTV <- as.numeric(inTaxa$PTV)

  samples <- unique(subset(inCts,select=c(sampID,'SAMPID')))

  inCts.1 <- inCts[inCts$TAXA_ID %in% inTaxa$TAXA_ID, c('SAMPID','TAXA_ID','TOTAL','IS_DISTINCT')]
  inCts.1 <- inCts.1[!is.na(inCts.1$TOTAL) & inCts.1$TOTAL>0,]

  inTaxa.1 <- inTaxa
  inTaxa.1$TOLR <- with(inTaxa.1, ifelse(PTV >= 7, 1, NA))
  inTaxa.1$FACL <- with(inTaxa.1, ifelse(PTV >3 & PTV < 7, 1, NA))
  inTaxa.1$INTL <- with(inTaxa.1, ifelse(PTV <= 3, 1, NA))
  inTaxa.1$NTOL <- with(inTaxa.1, ifelse(PTV < 6, 1, NA))
  inTaxa.1$STOL <- with(inTaxa.1, ifelse(PTV >= 8, 1, NA))
  inTaxa.1$TL01 <- with(inTaxa.1, ifelse(PTV < 2, 1, NA))
  inTaxa.1$TL23 <- with(inTaxa.1, ifelse(PTV >= 2 & PTV < 4, 1, NA))
  inTaxa.1$TL45 <- with(inTaxa.1, ifelse(PTV >= 4 & PTV < 6, 1, NA))
  inTaxa.1$TL67 <- with(inTaxa.1, ifelse(PTV >= 6 & PTV < 8, 1, NA))

  # Drop non-target taxa if included in taxalist
  if(length(grep('NON_TARGET',names(inTaxa.1)))>0) {
    inTaxa.1 <- subset(inTaxa.1,is.na(NON_TARGET)|NON_TARGET=='')
  }

  params<-c('TOLR','FACL','INTL','NTOL','STOL','TL01','TL23','TL45','TL67','PTV')

  taxalong <- reshape(inTaxa.1[,c('TAXA_ID',params)], idvar = 'TAXA_ID', direction = 'long',
                      varying = params, timevar = 'TRAIT', v.names = 'value',
                      times = params)

  taxalong <- taxalong[!is.na(taxalong$value),]

  taxalong$TRAIT <- as.character(taxalong$TRAIT)

  totals <- aggregate(x = list(TOTLNIND = inCts.1$TOTAL, TOTLNTAX = inCts.1$IS_DISTINCT), by = inCts.1[c('SAMPID')],
                      FUN = sum)

  inCts.1 <- merge(inCts.1, totals, by = 'SAMPID')
  inCts.1$CALCPIND <- with(inCts.1, TOTAL/TOTLNIND)
  inCts.1$CALCPTAX <- with(inCts.1, IS_DISTINCT/TOTLNTAX)

  # Merge the count data with the taxalist containing only the traits of
  # interest
  traitDF <- merge(inCts.1, taxalong, by='TAXA_ID')

  # Calculate no. individuals, % individuals, no. taxa, and % taxa for each
  # trait in taxalist
  outMet.1 <- aggregate(x = list(NIND = traitDF$TOTAL, NTAX = traitDF$IS_DISTINCT),
                        by = traitDF[c('SAMPID','TRAIT','TOTLNTAX')],
                        FUN = sum)
  outMet.2 <- aggregate(x = list(PIND = traitDF$CALCPIND, PTAX = traitDF$CALCPTAX),
                        by = traitDF[c('SAMPID','TRAIT','TOTLNTAX')],
                        FUN = function(x){round(sum(x)*100, 2)})

  outMet <- merge(outMet.1, outMet.2, by = c('SAMPID','TRAIT','TOTLNTAX'))

  # Melt df to create metric names, then recast into wide format with metric
  # names
  print("Done calculating tolerance metrics.")

  outLong <- reshape(outMet, idvar = c('SAMPID','TOTLNTAX','TRAIT'), direction = 'long',
                     varying = names(outMet)[!names(outMet) %in% c('SAMPID','TOTLNTAX','TRAIT')],
                     timevar = 'variable', v.names = 'value',
                     times = names(outMet)[!names(outMet) %in% c('SAMPID','TOTLNTAX','TRAIT')])

  outLong$variable <- paste(outLong$TRAIT,outLong$variable,sep='')
  outLong$TRAIT <- NULL

  outWide <- reshape(outLong, idvar = c('SAMPID','TOTLNTAX'), direction = 'wide',
                     timevar = 'variable', v.names = 'value')

  names(outWide) <- gsub("value\\.", "", names(outWide))

  outWide <- merge(outWide, samples, by='SAMPID')

  # Calculate HBI (WTD_TV)
  if(nrow(subset(taxalong,TRAIT %in% c('PTV')))>0){
    TVI <- tolindex(inCts.1,taxalong)
    outWide <- merge(outWide,TVI,by="SAMPID",all.x=TRUE)
    print("Calculated tolerance index.")
  }

  empty_tax <- data.frame(t(rep(NA,28)),stringsAsFactors=F)
  names(empty_tax) <- c('FACLNTAX','FACLPIND','FACLPTAX'
                        ,'INTLNTAX','INTLPIND','INTLPTAX'
                        ,'NTOLNTAX','NTOLPIND','NTOLPTAX'
                        ,'STOLNTAX','STOLPIND','STOLPTAX'
                        ,'TL01NTAX','TL01PIND','TL01PTAX'
                        ,'TL23NTAX','TL23PIND','TL23PTAX'
                        ,'TL45NTAX','TL45PIND','TL45PTAX'
                        ,'TL67NTAX','TL67PIND','TL67PTAX'
                        ,'TOLRNTAX','TOLRPIND','TOLRPTAX'
                        ,'WTD_TV')

  outWide.all <- merge(outWide, empty_tax, all=TRUE)
  outWide.all <- subset(outWide.all, !is.na(SAMPID))

  # # If we re-melt df now, we have missing values where the metric should be a
  # # zero, so we can set NAs to 0 now
  updNames <- names(outWide.all)[names(outWide.all) %nin% c('WTD_TV','NAT_WTD_TV','SAMPID')]
  outWide.all[,updNames] <- lapply(outWide.all[,updNames], function(x){ifelse(is.na(x), 0, x)})
  # # Finally, we can recast the metrics df into wide format for output
  outWide.fin <- outWide.all
  outWide.fin$SAMPID <- NULL

  outWide.fin <- outWide.fin[, c(sampID,'FACLNTAX','FACLPIND','FACLPTAX'
                                 ,'INTLNTAX','INTLPIND','INTLPTAX'
                                 ,'NTOLNTAX','NTOLPIND','NTOLPTAX'
                                 ,'STOLNTAX','STOLPIND','STOLPTAX'
                                 ,'TL01NTAX','TL01PIND','TL01PTAX'
                                 ,'TL23NTAX','TL23PIND','TL23PTAX'
                                 ,'TL45NTAX','TL45PIND','TL45PTAX'
                                 ,'TL67NTAX','TL67PIND','TL67PTAX'
                                 ,'TOLRNTAX','TOLRPIND','TOLRPTAX'
                                 ,'WTD_TV')]

  return(outWide.fin)

}


#' @export
#' @title Calculate dominance and diversity benthic metrics
#'
#' @description This function calculates overall dominance,
#' chironomid dominance, and Shannon diversity
#' benthic metrics.
#'
#' @param inCts A data frame containing, at minimum, the variables
#' specified in the arguments for sampID, dist, ct, and taxa_id
#' @param inTaxa a data frame containing taxonomic information,
#' including variables for PHYLUM, CLASS, ORDER, FAMILY, SUBFAMILY,
#' and TRIBE, as well as autecology traits with names that match those
#' in the arguments ffg, habit, and ptv. In addition, there
#' should be a variable with the name in argument taxa_id that matches
#' with all of those in the indf data frame
#' @param sampID A character vector containing the names of all
#' variables in indf that specify a unique sample. If not specified,
#' the default is \emph{UID}
#' @param dist A string with the name of the distinctness variable,
#' which is assumed to have only values of 0 or 1. If not specified,
#' the default is \emph{IS_DISTINCT}.
#' @param ct A string with the name of the count variable. If not
#' specified, the default is \emph{TOTAL}.
#' @param taxa_id A string with the name of the taxon ID variable
#' in \emph{inCts} that matches that in \emph{inTaxa}. The default
#' value is \emph{TAXA_ID}.
#' @return A data frame containing the variables in sampID and
#' the benthic macroinvertebrate metrics as additional variables.
#' The names of metrics include HPRIME, DOM1PIND, DOM3PIND,
#' DOM5PIND, CHIRDOM1PIND, CHIRDOM3PIND, and CHIRDOM5PIND.
#'
#' Metric descriptions are included in \emph{NRSA_Invertebrate_Metric_Descriptions.pdf},
#' included in this package.
#' @author Karen Blocksom \email{Blocksom.Karen@epa.gov}

calcBentDominMets <- function(inCts, inTaxa, sampID="UID", dist="IS_DISTINCT",
                          ct="TOTAL",taxa_id="TAXA_ID"){

  ctVars <- c(sampID,dist,ct,taxa_id)
  if(any(ctVars %nin% names(inCts))){
    msgTraits <- which(ctVars %nin% names(inCts))
    print(paste("Missing variables in input data frame:",paste(ctVars[msgTraits],collapse=',')))
    return(NULL)
  }

  inCts <- subset(inCts,select=c(sampID,ct,dist,taxa_id))
  # Rename ct and dist to TOTAL and IS_DISTINCT
  names(inCts)[names(inCts)==ct] <- 'TOTAL'
  names(inCts)[names(inCts)==dist] <- 'IS_DISTINCT'
  names(inCts)[names(inCts)==taxa_id] <- 'TAXA_ID'

  for(i in 1:length(sampID)){
    if(i==1) inCts$SAMPID <- inCts[,sampID[i]]
    else inCts$SAMPID <- paste(inCts$SAMPID,inCts[,sampID[i]],sep='.')
  }

  # Make sure all taxa match to taxalist and send error if not
  checkTaxa <- inCts[!(inCts$TAXA_ID %in% inTaxa$TAXA_ID),]
  if(nrow(checkTaxa)>0){
    return(print('Taxa in counts that do not have matches in taxalist! Cannot continue.'))
  }

  # If necessary, load the bentTaxa data frame and assign it to inTaxa.  Though
  # NON_TARGET taxa are included in the table provided by NRSA, we need to exclude
  # them from our calculations.
  if('NON_TARGET' %in% names(inTaxa)){
    inTaxa <- subset(inTaxa, is.na(NON_TARGET) | NON_TARGET == "" |NON_TARGET=='N')
  }

  if('FAMILY' %nin% names(inTaxa)){
    print("Missing variable ORDER from input taxa list. Will not calculate chironomid dominance metrics.")
  }

  inTaxa <- subset(inTaxa,select=names(inTaxa) %in% c('TAXA_ID','FAMILY'))

  samples <- unique(subset(inCts,select=c(sampID,'SAMPID')))
  inCts.1 <- inCts[inCts$TAXA_ID %in% inTaxa$TAXA_ID, c('SAMPID','TAXA_ID','TOTAL','IS_DISTINCT')]
  inCts.1 <- inCts.1[!is.na(inCts.1$TOTAL) & inCts.1$TOTAL>0,]

  # Calculate Shannon Diversity
  outMet <- ShanDiversity(inCts.1)
  print("Shannon Diversity calculated")

  # Calculate % dominant individuals in top N taxa
  outMet <- merge(outMet,Dominance(inCts.1,1),by="SAMPID",all.x=TRUE)
  outMet <- merge(outMet,Dominance(inCts.1,3),by="SAMPID",all.x=TRUE)
  outMet <- merge(outMet,Dominance(inCts.1,5),by="SAMPID",all.x=TRUE)

  # In cases where there are fewer taxa than the topN number, we need to fix NAs
  # by setting to 100%
  outMet$DOM3PIND <- with(outMet, ifelse(is.na(DOM3PIND) & !is.nan(DOM1PIND), 100, DOM3PIND))
  outMet$DOM5PIND <- with(outMet, ifelse(is.na(DOM5PIND) & !is.nan(DOM1PIND), 100, DOM5PIND))
  print("Dominance calculated")

  if('FAMILY' %in% names(inTaxa)){
    # Now calculate % dominant individuals in the top N taxa, but using totlnind
    # in sample as base of calculation
    chiroIn <- merge(inCts,inTaxa[,c('TAXA_ID','FAMILY')],by="TAXA_ID")
    chiroIn <- subset(chiroIn, FAMILY=='CHIRONOMIDAE', select=c('SAMPID','TAXA_ID','TOTAL','IS_DISTINCT'))
    chiroIn$CALC <- with(chiroIn, IS_DISTINCT*TOTAL)

    chiroTotal <- aggregate(x = list(TOTLDIST = chiroIn$CALC), by = chiroIn[c('SAMPID')],
                            FUN = sum)

    chiroIn.1 <- merge(chiroIn, chiroTotal, by=c('SAMPID'))

    chiro1 <- Dominance(chiroIn.1, topN=1)
    names(chiro1)[names(chiro1)=='DOM1PIND'] <- 'CHIRDOM1PIND'

    chiro3 <- Dominance(chiroIn.1, topN=3)
    names(chiro3)[names(chiro3)=='DOM3PIND'] <- 'CHIRDOM3PIND'

    chiro5 <- Dominance(chiroIn.1, topN=5)
    names(chiro5)[names(chiro5)=='DOM5PIND'] <- 'CHIRDOM5PIND'

    outMet <- merge(outMet, chiro1, by='SAMPID', all.x=TRUE)
    outMet <- merge(outMet, chiro3, by='SAMPID', all.x=TRUE)
    outMet <- merge(outMet, chiro5, by='SAMPID', all.x=TRUE)

    outMet$CHIRDOM1PIND <- with(outMet, ifelse(is.na(CHIRDOM1PIND),0,CHIRDOM1PIND))
    outMet$CHIRDOM3PIND <- with(outMet, ifelse(is.na(CHIRDOM3PIND) & CHIRDOM1PIND>0
                                               ,100,ifelse(is.na(CHIRDOM3PIND) & CHIRDOM1PIND==0,0,CHIRDOM3PIND)))
    outMet$CHIRDOM5PIND <- with(outMet, ifelse(is.na(CHIRDOM5PIND) & CHIRDOM3PIND>0, 100
                                               ,ifelse(is.na(CHIRDOM5PIND) & CHIRDOM3PIND==0, 0, CHIRDOM5PIND)))

    print("Chironomid dominance calculated")

  }

  empty_tax <- data.frame(t(rep(NA,7)),stringsAsFactors=F)
  names(empty_tax) <- c('HPRIME','DOM1PIND','DOM3PIND','DOM5PIND'
                        ,'CHIRDOM1PIND','CHIRDOM3PIND','CHIRDOM5PIND')

  outWide <- merge(outMet, empty_tax, all=TRUE)
  outWide <- subset(outWide, !is.na(SAMPID))
  outWide <- merge(outWide, samples, by='SAMPID')

  outWide[is.na(outWide)] <- 0

  outWide.fin <- outWide
  outWide.fin$SAMPID <- NULL

  outWide.fin <- outWide.fin[, c(sampID,'HPRIME','DOM1PIND','DOM3PIND','DOM5PIND'
                                 ,'CHIRDOM1PIND','CHIRDOM3PIND','CHIRDOM5PIND')]

  return(outWide.fin)
}
