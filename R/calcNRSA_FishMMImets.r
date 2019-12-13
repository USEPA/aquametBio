#' @export
#' @title Calculate only fish metrics necessary for NRSA MMI
#'
#' @description This is a function that calculates
#' only the fish metrics used in the National Rivers and Streams
#' Assessment fish MMI.
#'
#' @param indata A data frame containing, at minimum, the variables
#' specified in the arguments for \emph{sampID}, \emph{dist},
#' \emph{ct}, and \emph{taxa_id}, as well as the optional variable
#' for non-native status in \emph{nonnat}.
#' @param inTaxa a data frame containing taxonomic information,
#' including variables that match the arguments for (\emph{family}),
#' (\emph{genus}), and (\emph{comname}), as well as autecology traits
#' with names that match those in the arguments \emph{tol}, \emph{vel},
#' \emph{habitat}, \emph{trophic}, \emph{migr}, and \emph{reprod}.
#' In addition, there should be a variable with the
#' name in argument \emph{taxa_id} that matches
#' with all of those in the indata data frame
#' @param sampID A character vector containing the names of all
#' variables in indata that specify a unique sample. If not specified,
#' the default is \emph{UID}
#' @param ecoreg A string with the name of the ecoregion variable.
#' Valid values that correspond to regions used in NRSA are
#' CPL, NAP, NPL, SAP, sPL, TPL, UMW, WMT, and XER.
#' @param dist A string with the name of the distinctness variable,
#' which is assumed to have only values of 0 or 1. If not specified,
#' the default is \emph{IS_DISTINCT}.
#' @param ct A string with the name of the count variable. If not
#' specified, the default is \emph{TOTAL}.
#' @param taxa_id A string with the name of the taxon ID variable
#' in \emph{indata} that matches that in \emph{inTaxa}. The default
#' value is \emph{TAXA_ID}.
#' @param tol A string with the name of the variable in the
#' \emph{inTaxa} taxalist containing tolerance categories. Valid
#' values include S (sensitive), I (intermediate), and T (tolerant).
#' The default value is \emph{TOLERANCE}.
#' @param vel A strings with the name of the variable in the
#' \emph{inTaxa} taxalist containing velocity preference values.
#' Valid values include R (Rheophil), P (Pool), O (Other), or
#' blank, if unknown. The default value is \emph{VELOCITY}.
#' @param habitat A string with the name of the variable in
#' \emph{inTaxa} containing habitat preference values. Valid
#' values include B (benthic), W (water column), E (edge), or
#' blank, if unknown. The default value is \emph{HABITAT}.
#' @param trophic A string with the name of the variable in
#' \emph{inTaxa} containing trophic category values. Valid
#' values include I (invertivore), C (carnivore), O (omnivore),
#' H (herbivore), or blank if unknown. The default value is
#' \emph{inTaxa} is TROPHIC.
#' @param migr A string with the name of the variable in
#' \emph{inTaxa} containing migratory status value. Valid
#' values include N (No), Y (Yes), and blank if unknown. The
#' default value is \emph{MIGRATORY}.
#' @param nonnat A string with the name of the optional variable in
#' \emph{inCts} containing non-native status. Valid values are 'Y' for
#' non-native and 'N' for native. The default name
#' is \emph{NONNATIVE}.
#' @param reprod A string with the name of the variable in
#' \emph{inTaxa} containing reproductive trait values. Valid
#' values include C (Clean, coarse lithophil), D (Drifter),
#' G (Guarder), O (Other), or blank if unknown. The default
#' value is \emph{REPROD}.
#' @param family A string with the name of the variable in the
#' \emph{inTaxa} taxalist containing family name. The default value
#' is \emph{FAMILY}.
#' @param genus A string with the name of the variable in the
#' \emph{inTaxa} taxalist containing the genus name. The default
#' value is \emph{GENUS}
#' @param comname A string with the name of the variable in the
#' \emph{inTaxa} taxalist containing the common name. The
#' default value is \emph{NAME}.
#' @return A data frame containing the variables in sampID,
#' the total number of individuals in the sample as TOTLNIND, and
#' all fish metrics used in MMIs as additional variables. All metrics
#' will be separate columns, but values will only be provided for
#' the metrics in the MMI for the specific region of each site,
#' with all other metrics being missing for that site. The metrics,
#' by aggregated ecoregion, are:
#'
#'  CPL: ALIENPIND, RBCATONTAX, LOTPIND, INTLMIGRPTAX, LITHPIND, NAT_TOTLNTAX,
#'      TOLRNTAX, INVPTAX
#'
#'  NAP: ALIENNTAX, SALMNTAX, NAT_RHEOPIND, INTLMIGRPIND, LITHPTAX,
#'      NTOLPTAX, TOLRNTAX, INVNTAX
#'
#'  NPL: ALIENNTAX, NAT_CYPRPIND, LOTNTAX, MIGRNTAX, LITHPIND, NTOLPTAX,
#'      NAT_INTLPIND, NAT_CARNNTAX
#'
#'  SAP: NAT_PTAX, NAT_CENTNTAX, NAT_NTOLBENTPTAX, NAT_MIGRNTAX,
#'      NAT_LITHPIND, NTOLPTAX, TOLRPTAX, INVPIND
#'
#'  SPL: NAT_PIND, CYPRPTAX, RHEOPIND, NAT_MIGRPTAX, LITHNTAX,
#'      NAT_NTOLNTAX, TOLRNTAX, HERBPTAX
#'
#'  TPL: ALIENNTAX, NAT_ICTAPIND, RHEONTAX, INTLMIGRNTAX, LITHPIND,
#'      NAT_NTOLNTAX, INTLPTAX, CARNNTAX
#'
#'  UMW: NAT_PTAX, CYPRNTAX, INTLLOTNTAX, INTLMIGRPTAX, LITHPIND,
#'      NTOLNTAX, TOLRNTAX, INTLINVPTAX
#'
#'  WMT: NAT_PIND, NAT_CATOPIND, INTLLOTPTAX, NAT_MIGRPTAX, LITHPTAX,
#'      NAT_TOTLNTAX, TOLRNTAX, NAT_HERBPTAX
#'
#'  XER: NAT_PIND, CENTPTAX, RHEOPIND, MIGRPTAX, LITHNTAX, NTOLPTAX,
#'      TOLRNTAX, BENTINVPTAX
#'
#' Metric descriptions are included in \emph{NRSA_Fish_Metric_Descriptions.pdf},
#' included in this package.
#' @author Karen Blocksom \email{Blocksom.Karen@epa.gov}
#' @keywords survey
calcNRSA_FishMMImets <- function(indata,inTaxa=NULL, sampID="UID", ecoreg=NULL
                                 ,dist="IS_DISTINCT",ct="TOTAL"
                                 ,taxa_id='TAXA_ID',tol='TOLERANCE'
                                 ,vel='VELOCITY', habitat='HABITAT'
                                 ,trophic='TROPHIC', migr='MIGRATORY', nonnat='NONNATIVE'
                                 ,reprod='REPROD', family='FAMILY', genus='GENUS'
                                 ,comname='NAME'){

  if(is.null(inTaxa)) {
    inTaxa <- fishTaxa
    inTaxa <- subset(inTaxa, is.na(NON_TARGET) | NON_TARGET == "")
  }

  ctVars <- c(sampID, dist, ct, taxa_id, nonnat,ecoreg)
  if(any(ctVars %nin% names(indata))){
    msgTraits <- which(ctVars %nin% names(indata))
    print(paste("Missing variables in input count data frame:",paste(ctVars[msgTraits],collapse=',')))
    return(NULL)
  }

  ecoCk <- unique(indata[,ecoreg])
  ecos <- c('CPL','NAP','NPL','SAP','SPL','TPL','UMW','WMT','XER')
  if(any(ecoCk %nin% ecos)){
    msgEco <- which(ecoCk %nin% ecos)
    print(paste("These ecoregions are not valid: "
                ,paste(ecoCk[msgEco],collapse=',')))
    return(NULL)
  }

  # Combine all values in sampID into one sampID in df
  for(i in 1:length(sampID)){
    if(i==1) indata$SAMPID <- indata[,sampID[i]]
    else indata$SAMPID <- paste(indata$SAMPID,indata[,sampID[i]],sep='.')
  }
  # Keep data frame with crosswalk info between sampID and SAMPID
  samples <- unique(subset(indata,select=c(sampID,ecoreg,'SAMPID')))


  metnames <- data.frame(ECO=c(rep('CPL',8),rep('NAP',8),rep('NPL',8),rep('SAP',8),rep('SPL',8)
                                          ,rep('TPL',8),rep('UMW',8),rep('WMT',8),rep('XER',8))
                         ,METRIC=c(c('ALIENPIND','RBCATONTAX','LOTPIND','INTLMIGRPTAX','LITHPIND','NAT_TOTLNTAX','TOLRNTAX'
                                        ,'INVPTAX')
                                      ,c('ALIENNTAX','SALMNTAX','NAT_RHEOPIND','INTLMIGRPIND','LITHPTAX','NTOLPTAX','TOLRNTAX'
                                         ,'INVNTAX')
                                      ,c('ALIENNTAX','NAT_CYPRPIND','LOTNTAX','MIGRNTAX','LITHPIND','NTOLPTAX','NAT_INTLPIND'
                                         ,'NAT_CARNNTAX')
                                      ,c('NAT_PTAX','NAT_CENTNTAX','NAT_NTOLBENTPTAX','NAT_MIGRNTAX','NAT_LITHPIND','NTOLPTAX'
                                         ,'TOLRPTAX','INVPIND')
                                      ,c('NAT_PIND','CYPRPTAX','RHEOPIND','NAT_MIGRPTAX','LITHNTAX','NAT_NTOLNTAX','TOLRNTAX'
                                         ,'HERBPTAX')
                                      ,c('ALIENNTAX','NAT_ICTAPIND','RHEONTAX','INTLMIGRNTAX','LITHPIND','NAT_NTOLNTAX','INTLPTAX'
                                         ,'CARNNTAX')
                                      ,c('NAT_PTAX','CYPRNTAX','INTLLOTNTAX','INTLMIGRPTAX','LITHPIND','NTOLNTAX','TOLRNTAX'
                                         ,'INTLINVPTAX')
                                      ,c('NAT_PIND','NAT_CATOPIND','INTLLOTPTAX','NAT_MIGRPTAX','LITHPTAX','NAT_TOTLNTAX'
                                         ,'TOLRNTAX','NAT_HERBPTAX')
                                      ,c('NAT_PIND','CENTPTAX','RHEOPIND','MIGRPTAX','LITHNTAX','NTOLPTAX','TOLRNTAX'
                                         ,'BENTINVPTAX'))
                         ,stringsAsFactors=F)

  empty_base <- data.frame(t(rep(NA,length(unique(metnames$METRIC)))),stringsAsFactors=F)
  names(empty_base) <- unique(metnames$METRIC)

  # Taxonomy and traits checks
  necTraits <- c(taxa_id,tol, vel, habitat, trophic, migr,
                 reprod, family, genus, comname)
  if(any(necTraits %nin% names(inTaxa))){
    msgTraits <- which(necTraits %nin% names(inTaxa))
    print(paste("Some of the traits are missing from the taxa list. The following are required for metric calculations to run:"
                , necTraits[msgTraits]))
    return(NULL)
  }

  inTaxa <- subset(inTaxa,select=names(inTaxa) %in% necTraits)

  # Rename counts and distinct variables to TOTAL and IS_DISTINCT
  names(indata)[names(indata)==ct] <- 'TOTAL'
  names(indata)[names(indata)==dist] <- 'IS_DISTINCT'
  names(indata)[names(indata)==taxa_id] <- 'TAXA_ID'
  names(indata)[names(indata)==nonnat] <- 'NONNATIVE'
  names(inTaxa)[names(inTaxa)==taxa_id] <- 'TAXA_ID'

  names(inTaxa)[names(inTaxa)==tol] <- 'TOLERANCE'
  names(inTaxa)[names(inTaxa)==vel] <- 'VELOCITY'
  names(inTaxa)[names(inTaxa)==habitat] <- 'HABITAT'
  names(inTaxa)[names(inTaxa)==trophic] <- 'TROPHIC'
  names(inTaxa)[names(inTaxa)==migr] <- 'MIGRATORY'
  names(inTaxa)[names(inTaxa)==family] <- 'FAMILY'
  names(inTaxa)[names(inTaxa)==genus] <- 'GENUS'
  names(inTaxa)[names(inTaxa)==comname] <- 'NAME'
  names(inTaxa)[names(inTaxa)==reprod] <- 'REPROD'

  indata[,c('TOTAL','IS_DISTINCT')] <- lapply(indata[,c('TOTAL','IS_DISTINCT')],as.numeric)
  indata$TAXA_ID <- as.character(indata$TAXA_ID)
  inTaxa$TAXA_ID <- as.character(inTaxa$TAXA_ID)

  ## for inCts1, keep only observations without missing or zero TOTAL values or TAXA_ID and TAXA_ID!=99999
  indata.1 <- subset(indata,!is.na(TAXA_ID) & !is.na(TOTAL) & TOTAL!=0)

  indata.2 <- subset(indata.1, select = c('SAMPID','TAXA_ID','TOTAL','IS_DISTINCT','NONNATIVE'))

  inTaxa.1 <- inTaxa
  inTaxa.1$BENTINV <- with(inTaxa.1, ifelse(HABITAT=='B' & TROPHIC=='I',1,NA))
  inTaxa.1$CARN <- with(inTaxa.1, ifelse(TROPHIC=='C',1,NA))
  inTaxa.1$CARN <- with(inTaxa.1, ifelse(TROPHIC=='C',1,NA))
  inTaxa.1$CENT <- with(inTaxa.1, ifelse(toupper(FAMILY)=='CENTRARCHIDAE' & toupper(GENUS)!='MICROPTERUS',1,NA))
  inTaxa.1$CYPR <- with(inTaxa.1, ifelse(toupper(FAMILY)=='CYPRINIDAE' & toupper(NAME) %nin%
                 c('COMMON CARP','GOLDFISH','BIGHEAD CARP','GRASS CARP','MIRROR CARP'),1,NA))
  inTaxa.1$HERB <- with(inTaxa.1, ifelse(TROPHIC=='H',1,NA))
  inTaxa.1$INTL <- with(inTaxa.1, ifelse(TOLERANCE=='S',1,NA))
  inTaxa.1$INTLINV <- with(inTaxa.1, ifelse(INTL==1 & TROPHIC=='I',1,NA))
  inTaxa.1$INTLLOT <- with(inTaxa.1, ifelse(INTL==1 & VELOCITY %in% c('R','O'),1,NA))
  inTaxa.1$INTLMIGR <- with(inTaxa.1, ifelse(INTL==1 & MIGRATORY=='Y',1,NA))
  inTaxa.1$INV <- with(inTaxa.1, ifelse(TROPHIC %in% c('I'),1,NA))
  inTaxa.1$LITH <- with(inTaxa.1, ifelse(REPROD=='C',1,NA))
  inTaxa.1$LOT <- with(inTaxa.1, ifelse(VELOCITY %in% c('R','O'),1,NA))
  inTaxa.1$MIGR <- with(inTaxa.1, ifelse(MIGRATORY=='Y',1,NA))
  inTaxa.1$CATO <- with(inTaxa.1, ifelse(toupper(FAMILY)=='CATOSTOMIDAE',1,NA))
  inTaxa.1$ICTA <- with(inTaxa.1, ifelse(toupper(FAMILY)=='ICTALURIDAE',1,NA))
  inTaxa.1$NTOL <- with(inTaxa.1, ifelse(TOLERANCE %in% c('S','I'),1,NA))
  inTaxa.1$NTOLBENT <- with(inTaxa.1, ifelse(NTOL==1 & HABITAT=='B',1,NA))
  inTaxa.1$RHEO <- with(inTaxa.1, ifelse(VELOCITY=='R',1,NA))
  inTaxa.1$RBCATO <- with(inTaxa.1, ifelse(toupper(GENUS) %in% c('MOXOSTOMA', 'HYPENTELIUM'
                                       , 'MINYTREMA', 'ERIMYZON' , 'CATOSTOMUS', 'CYCLEPTUS'
                                       , 'PANTOSTEUS' , 'THOBURNIA'),1,NA))
  inTaxa.1$SALM <- with(inTaxa.1, ifelse(toupper(FAMILY)=='SALMONIDAE',1,NA))
  inTaxa.1$TOLR <- with(inTaxa.1, ifelse(TOLERANCE=='T',1,NA))

  params<-c('BENTINV','CARN','CENT','CYPR','HERB','INV','INTLINV','INTLLOT','INTLMIGR','INTL','LITH','LOT','MIGR','CATO'
          ,'ICTA','NTOLBENT','NTOL','RHEO','RBCATO','SALM','TOLR')

  inTaxa.2 <- subset(inTaxa.1,select=names(inTaxa.1) %in% c('TAXA_ID',params))

  params.use <- names(inTaxa.2)[names(inTaxa.2) %in% params]

  taxalong <- reshape(inTaxa.2, idvar = c('TAXA_ID'), direction = 'long',
                      varying = params.use, times = params.use, v.names = 'value', timevar = 'TRAIT')
  taxalong <- subset(taxalong, !is.na(value))

  totals <- aggregate(x = list(TOTLNIND = indata.2$TOTAL, TOTLNTAX = indata.2$IS_DISTINCT), by = indata.2[c('SAMPID')],
                      FUN = sum)

  indata.3 <- merge(indata.2, totals, by = 'SAMPID')
  indata.3$CALCPIND <- with(indata.3, TOTAL/TOTLNIND)
  indata.3$CALCPTAX <- with(indata.3, IS_DISTINCT/TOTLNTAX)

  # Merge the count data with the taxalist containing only the traits of
  # interest
  traitDF <- merge(indata.3, taxalong, by='TAXA_ID')

  # Calculate no. individuals, % individuals, no. taxa, and % taxa for each
  # trait in taxalist
  outMet.1 <- aggregate(x = list(NTAX = traitDF$IS_DISTINCT),
                        by = traitDF[c('SAMPID','TRAIT')],
                        FUN = sum)
  outMet.2 <- aggregate(x = list(PIND = traitDF$CALCPIND, PTAX = traitDF$CALCPTAX),
                        by = traitDF[c('SAMPID','TRAIT')],
                        FUN = function(x){round(sum(x)*100, 2)})

  outMet <- merge(outMet.1, outMet.2, by = c('SAMPID','TRAIT'))

  # Melt df to create metric names, then recast into wide format with metric
  # names
  outLong <- reshape(outMet, idvar = c('SAMPID','TRAIT'), direction = 'long',
                     varying = c('PIND','PTAX','NTAX'), timevar = 'variable',
                     v.names = 'value', times = c('PIND','PTAX','NTAX'))

  outLong$variable <- paste(outLong$TRAIT,outLong$variable,sep='')
  outLong$TRAIT <- NULL

  outWide <- reshape(outLong, idvar = 'SAMPID', direction = 'wide',
                     v.names = 'value', timevar = 'variable')
  names(outWide) <- gsub("value\\.", "", names(outWide))

  outWide <- merge(outWide, totals, by = 'SAMPID', all.y=TRUE)

  # Calculate metrics based on native status and other traits
  if(any(unique(indata.3$NON_NATIVE) %nin% c('Y','N'))){
    return(print("No native and alien datasets were created because NON_NATIVE must only be 'Y' or 'N' values"))
  }else{
    inNative <- subset(indata.3,NONNATIVE=='N')
    if(length(inNative)>0){
      natTot <- aggregate(x = list(NAT_TOTLNIND = inNative$TOTAL, NAT_TOTLNTAX = inNative$IS_DISTINCT),
                          by = inNative[c('SAMPID')], FUN = sum)

      inNative.tot <- merge(inNative, natTot, by = 'SAMPID')

      natMets <- merge(inNative.tot, taxalong, by = 'TAXA_ID')
      natMets$CALCPIND <- with(natMets, TOTAL/NAT_TOTLNIND)
      natMets$CALCPTAX <- with(natMets, IS_DISTINCT/NAT_TOTLNTAX)

      natMets.1 <- aggregate(x = list(NTAX = natMets$IS_DISTINCT),
                             by = natMets[c('SAMPID','TRAIT')],
                             FUN = sum)
      natMets.2 <- aggregate(x = list(PIND = natMets$CALCPIND, PTAX = natMets$CALCPTAX),
                             by = natMets[c('SAMPID','TRAIT')],
                             FUN = function(x){round(sum(x)*100, 2)})

      natMets.comb <- merge(natMets.1, natMets.2, by = c('SAMPID','TRAIT'))

      natMets.long <- reshape(natMets.comb, idvar = c('SAMPID','TRAIT'), direction = 'long',
                              varying = c('NTAX', 'PIND', 'PTAX'), timevar = 'variable', v.names = 'value',
                              times = c('NTAX', 'PIND', 'PTAX'))
      natMets.long$variable <- with(natMets.long, paste('NAT_',TRAIT,variable,sep=''))
      natMets.long$TRAIT <- NULL

      natMets.fin <- reshape(natMets.long, idvar = c('SAMPID'), direction = 'wide',
                             v.names = 'value', timevar = 'variable')
      names(natMets.fin) <- gsub("value\\.", "", names(natMets.fin))

      outWide.1 <- merge(outWide, natMets.fin, all=TRUE)

    }else{
      outWide.1 <- outWide
    }
  }


  # Metrics relying only on native status

  inNat <- indata.3
  inNat$ALIEN <- with(inNat, ifelse(NONNATIVE=='Y',1,NA))
  inNat$NAT <- with(inNat, ifelse(NONNATIVE=='N',1,NA))
  inNat$ALIENNTAX <- with(inNat, IS_DISTINCT*ALIEN)
  inNat$ALIENPIND <- with(inNat, TOTAL*ALIEN/TOTLNIND)
  inNat$NAT_NTAX <- with(inNat, IS_DISTINCT*NAT)
  inNat$NAT_PIND <- with(inNat, TOTAL*NAT/TOTLNIND)
  inNat$NAT_PTAX <- with(inNat, IS_DISTINCT*NAT/TOTLNTAX)

  outNat.1 <- aggregate(x = list(ALIENNTAX=inNat$ALIENNTAX, NAT_TOTLNTAX=inNat$NAT_NTAX), by = inNat[c('SAMPID')],
                   FUN = function(x){sum(x, na.rm=T)})

  outNat.2 <- aggregate(x = list(ALIENPIND = inNat$ALIENPIND, NAT_PIND = inNat$NAT_PIND, NAT_PTAX = inNat$NAT_PTAX),
                   by = inNat[c('SAMPID')], FUN = function(x){round(sum(x, na.rm = TRUE)*100, 2)})

  outNat <- merge(outNat.1, outNat.2, by = 'SAMPID')

  outAll <- merge(outWide.1, outNat, by = 'SAMPID', all=TRUE)
  outAll <- merge(outAll, empty_base, all=TRUE)
  outAll <- subset(outAll, !is.na(SAMPID))
  outAll <- merge(outAll, samples, by='SAMPID', all.y=TRUE)

  varLong <- names(outAll)[names(outAll) %nin% c(sampID, ecoreg, 'SAMPID')]
  outLong <- reshape(outAll, idvar = c(sampID, ecoreg, 'SAMPID'), direction = 'long',
                     varying = varLong, times = varLong, timevar = 'variable', v.names = 'value')
  outLong$value <- with(outLong, ifelse(is.na(value),0,value))
  outLong <- merge(outLong, metnames, by.x = c(ecoreg, 'variable'), by.y = c('ECO', 'METRIC'))

  ckMetnum <- as.data.frame(table(SAMPID=outLong$SAMPID))
  ckMetnum <- subset(ckMetnum, Freq!=8)
  if(nrow(ckMetnum)>0){
    print("Error in output! Wrong number of metrics per site!")
    print(ckMetnum)
  }

  # Finally, we can recast the metrics df into wide format for output
  metOut <- reshape(outLong, idvar = c(sampID, 'SAMPID', ecoreg), direction = 'wide',
                    v.names = 'value', timevar = 'variable')
  names(metOut) <- gsub("value\\.", "", names(metOut))

  metOut <- merge(metOut, unique(indata.3[,c('SAMPID','TOTLNIND')]), by='SAMPID')
  metOut$SAMPID <- NULL

  return(metOut)
}




