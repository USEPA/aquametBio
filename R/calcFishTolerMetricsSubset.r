#' @export
#' @title Calculate all tolerance-related metrics
#' @description This function calculates all of the tolerance metrics,
#' and if additional trait values are included in the taxalist (habitat,
#' trophic, migr).
#' @param indata Input data frame containing variables as identified
#' in the arguments for \emph{sampID}, \emph{dist}, \emph{ct},
#' \emph{taxa_id}. If a variable for anomaly counts, as specified in
#' the argument \emph{anom}, is included, additional metrics are
#' calculated. The variable for the optional argument \emph{nonnat}
#' should also be included in this data frame.
#' @param inTaxa Data frame containing fish taxalist, along with autecology
#' traits. At a minimum, this taxalist must contain variables matching
#' the argument for \emph{tol}. If additional traits are included, as
#' specified by the arguments \emph{habitat}, \emph{trophic}, \emph{migr},
#' additional metrics that combine tolerance and other traits are also
#' calculated. The variable specified in the argument \emph{taxa_id}
#' must match that in \emph{indata}.
#' @param sampID sampID A character vector containing the names of all
#' variables in indf that specify a unique sample. If not specified,
#' the default is \emph{UID}.
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
#' @param tolval A string with the name of the variable in the
#' \emph{inTaxa} taxalist containing numeric tolerance values.
#' Valid values range from 0 to 10. The default value is TOL_VAL.
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
#' @return A data frame containing the variables in sampID and
#' the fish tolerance metrics as additional variables. Metric
#' descriptions are included in \emph{NRSA_Fish_Metric_Descriptions.pdf},
#' included in this package. The names of
#' metrics include INTLNIND, INTLNTAX, INTLPIND, INTLPTAX, MTOLNIND,
#' MTOLNTAX, MTOLPIND, MTOLPTAX, NTOLNIND, NTOLNTAX, NTOLPIND,
#' NTOLPTAX, TOLRNIND, TOLRNTAX, TOLRPIND, TOLRPTAX, WTD_TV, TOTLNIND,
#' and TOTLNTAX.
#'
#' If a non-native status variable is included, these metrics are also
#' calculated:
#' NAT_INTLNIND, NAT_INTLNTAX, NAT_INTLPIND, NAT_INTLPTAX, NAT_MTOLNIND,
#' NAT_MTOLNTAX, NAT_MTOLPIND, NAT_MTOLPTAX, NAT_NTOLNIND, NAT_NTOLNTAX,
#' NAT_NTOLPIND, NAT_NTOLPTAX,  NAT_TOLRNIND, NAT_TOLRNTAX, NAT_TOLRPIND,
#' NAT_TOLRPTAX, NAT_WTD_TV, NAT_TOTLNTAX, NAT_TOTLNIND,
#' NAT_PIND, NAT_PTAX.
#'
#' Additional metrics calculated if the appropriate additional traits
#' are included in \emph{inTaxa}: INTLINVNIND, INTLINVNTAX, INTLINVPIND,
#' INTLINVPTAX, INTLLOTNIND, INTLLOTNTAX, INTLLOTPIND, INTLLOTPTAX,
#' INTLMIGRNIND, INTLMIGRNTAX, INTLMIGRPIND, INTLMIGRPTAX,  INTLRHEONIND,
#' INTLRHEONTAX, INTLRHEOPIND, INTLRHEOPTAX, NTOLBENTNIND, NTOLBENTNTAX,
#' NTOLBENTPIND, NTOLBENTPTAX, NTOLCARNNIND, NTOLCARNNTAX,
#' NTOLCARNPIND, NTOLCARNPTAX, NTOLINVNIND, NTOLINVNTAX, NTOLINVPIND,
#' NTOLINVPTAX.
#'
#' If a non-native status variable is included, these additional metrics
#' are also calculated: NAT_INTLINVNIND, NAT_INTLINVNTAX, NAT_INTLINVPIND,
#' NAT_INTLINVPTAX, NAT_INTLLOTNIND, NAT_INTLLOTNTAX, NAT_INTLLOTPIND,
#' NAT_INTLLOTPTAX, NAT_INTLMIGRNIND, NAT_INTLMIGRNTAX, NAT_INTLMIGRPIND,
#' NAT_INTLMIGRPTAX, NAT_INTLRHEONIND, NAT_INTLRHEONTAX, NAT_INTLRHEOPIND,
#' NAT_INTLRHEOPTAX, NAT_NTOLBENTNIND, NAT_NTOLBENTNTAX, NAT_NTOLBENTPIND,
#' NAT_NTOLBENTPTAX, NAT_NTOLCARNNIND, NAT_NTOLCARNNTAX, NAT_NTOLCARNPIND,
#' NAT_NTOLCARNPTAX, NAT_NTOLINVNIND, NAT_NTOLINVNTAX, NAT_NTOLINVPIND,
#' NAT_NTOLINVPTAX.
#'
#' @author Karen Blocksom \email{Blocksom.Karen@epa.gov}
#'
calcFishTolMets <- function(indata, inTaxa=NULL, sampID='UID', dist='IS_DISTINCT'
                            , ct='TOTAL', taxa_id='TAXA_ID', tol='TOLERANCE'
                            , tolval='TOL_VAL', vel='VELOCITY', habitat='HABITAT'
                            , trophic='TROPHIC', migr='MIGRATORY', nonnat='NONNATIVE'){

  ctVars <- c(sampID,dist,ct,taxa_id)
  if(any(ctVars %nin% names(indata))){
    msgTraits <- which(ctVars %nin% names(indata))
    print(paste("Missing variables in input data frame:",paste(ctVars[msgTraits],collapse=',')))
    return(NULL)
  }

  if(nonnat %nin% names(indata)){
    print(paste("Cannot calculate native status-based metrics without a non-native status variable. Will calculate other metrics."))
  }

  # Combine all values in sampID into one sampID in df
  for(i in 1:length(sampID)){
    if(i==1) indata$SAMPID <- indata[,sampID[i]]
    else indata$SAMPID <- paste(indata$SAMPID,indata[,sampID[i]],sep='.')
  }
  # Keep data frame with crosswalk info between sampID and SAMPID
  samples <- unique(subset(indata,select=c(sampID,'SAMPID')))

  # If inTaxa is not specified, the default is the included fishTaxa dataset
  if(is.null(inTaxa)) {
    inTaxa <- fishTaxa
  }

  # Taxonomy and traits checks
  necTraits <- c(tol,tolval)
  if(any(necTraits %nin% names(inTaxa))){
    msgTraits <- which(necTraits %nin% names(inTaxa))
    return(paste("Some of the traits are missing from the taxa list. The following are required for metric calculations to run:", necTraits[msgTraits]))
  }
  optTraits <- c(habitat,trophic,migr,vel)
  if(any(optTraits %nin% names(inTaxa))){
    msgTraits <- which(optTraits %nin% names(inTaxa))
    print(paste("Optional traits are missing from the taxa list. Any tolerance metrics also using these traits will not be calculated:",
                 paste(optTraits[msgTraits],collapse=',')))
  }


  inTaxa <- subset(inTaxa,select=names(inTaxa) %in% c(taxa_id,tol,tolval,vel,habitat,trophic,migr))

  # Rename counts and distinct variables to TOTAL and IS_DISTINCT
  names(indata)[names(indata)==ct] <- 'TOTAL'
  names(indata)[names(indata)==dist] <- 'IS_DISTINCT'
  names(indata)[names(indata)==taxa_id] <- 'TAXA_ID'
  names(indata)[names(indata)==nonnat] <- 'NONNATIVE'
  names(inTaxa)[names(inTaxa)==taxa_id] <- 'TAXA_ID'

  names(inTaxa)[names(inTaxa)==tol] <- 'TOLERANCE'
  names(inTaxa)[names(inTaxa)==tolval] <- 'TOL_VAL'

  if(vel %in% names(inTaxa)){
    names(inTaxa)[names(inTaxa)==vel] <- 'VELOCITY'
  }

  if(habitat %in% names(inTaxa)){
    names(inTaxa)[names(inTaxa)==habitat] <- 'HABITAT'
  }

  if(trophic %in% names(inTaxa)){
    names(inTaxa)[names(inTaxa)==trophic] <- 'TROPHIC'
  }

  if(migr %in% names(inTaxa)){
    names(inTaxa)[names(inTaxa)==migr] <- 'MIGRATORY'
  }

  indata[,c('TOTAL','IS_DISTINCT')] <- lapply(indata[,c('TOTAL','IS_DISTINCT')],as.numeric)
  indata$TAXA_ID <- as.character(indata$TAXA_ID)
  inTaxa$TAXA_ID <- as.character(inTaxa$TAXA_ID)
  inTaxa$TOL_VAL <- as.numeric(inTaxa$TOL_VAL)

  ## for inCts1, keep only observations without missing or zero TOTAL values or TAXA_ID and TAXA_ID!=99999
  indata.1 <- subset(indata,!is.na(TAXA_ID) & !is.na(TOTAL) & TOTAL!=0)

  ## Now sum by TAXA_ID for ANOM_CT and for TOTAL for each sample
  # Two approaches depending on whether or not NON_NATIVE occurs in the counts data frame
  if('NONNATIVE' %in% names(indata.1)){
    maxDist <- aggregate(x = list(IS_DISTINCT = indata.1$IS_DISTINCT), by = indata.1[c('SAMPID','TAXA_ID','NONNATIVE')],
                         FUN = function(x){max(as.integer(x))})
    sumTot <- aggregate(x = list(TOTAL = indata.1$TOTAL), by = indata.1[c('SAMPID','TAXA_ID','NONNATIVE')], FUN = sum)

    indata.2 <- merge(maxDist, sumTot, by = c('SAMPID','TAXA_ID','NONNATIVE'))

    # indata.2 <- plyr::ddply(indata.1,c('SAMPID','TAXA_ID','NONNATIVE'),summarise,IS_DISTINCT=max(IS_DISTINCT),TOTAL=sum(TOTAL))
    CALCNAT <- 'Y'
  }else{
    maxDist <- aggregate(x = list(IS_DISTINCT = indata.1$IS_DISTINCT), by = indata.1[c('SAMPID','TAXA_ID')],
                         FUN = function(x){max(as.integer(x))})
    sumTot <- aggregate(x = list(TOTAL = indata.1$TOTAL), by = indata.1[c('SAMPID','TAXA_ID')], FUN = sum)

    indata.2 <- merge(maxDist, sumTot, by = c('SAMPID','TAXA_ID'))

        # indata.2 <- plyr::ddply(indata.1,c('SAMPID','TAXA_ID'),summarise,IS_DISTINCT=max(IS_DISTINCT),TOTAL=sum(TOTAL))
    CALCNAT <- 'N'
  }
  # Find all samples with a missing TAXA_ID, which means there are no counts for the site, and output the rows so the user can
  ## verify no sample was collected
  if(nrow(subset(indata.2,is.na(TAXA_ID)))>0) {
    print("Make sure these missing TAXA_ID values are valid.")
    print(subset(indata,is.na(TAXA_ID)))
  }

  # Make sure all necessary columns in inCts are numeric
  inCts <- indata.2
  inCts$TOTAL <- with(inCts, as.numeric(TOTAL))
  inCts$IS_DISTINCT <- with(inCts, as.integer(IS_DISTINCT))
  # inCts <- plyr::mutate(indata.2,TOTAL=as.numeric(TOTAL),IS_DISTINCT=as.integer(IS_DISTINCT))

  inCts.1 <- inCts[inCts$TAXA_ID %in% inTaxa$TAXA_ID,]
  # inCts.1 <- dplyr::semi_join(inCts,subset(inTaxa,select='TAXA_ID'),by='TAXA_ID')

  if(CALCNAT=='Y'){
    inCts.1 <- inCts.1[!is.na(inCts.1$TOTAL) & inCts.1$TOTAL>0,c('SAMPID','TAXA_ID','TOTAL','IS_DISTINCT','NONNATIVE')]
    # inCts.1 <- dplyr::select(inCts.1,SAMPID, TAXA_ID, TOTAL, IS_DISTINCT,NONNATIVE) %>%
    # subset(!is.na(TOTAL) & TOTAL>0)
  }else{
    inCts.1 <- inCts.1[!is.na(inCts.1$TOTAL) & inCts.1$TOTAL>0,c('SAMPID','TAXA_ID','TOTAL','IS_DISTINCT')]
    # inCts.1 <- dplyr::select(inCts.1,SAMPID, TAXA_ID, TOTAL, IS_DISTINCT) %>%
    #   subset(!is.na(TOTAL) & TOTAL>0)
  }
  # Now create indicator variables
  inTaxa.1 <- inTaxa
  inTaxa.1$NTOL <- with(inTaxa.1, ifelse(TOLERANCE %in% c('S','I'),1,NA))
  inTaxa.1$INTL <- with(inTaxa.1, ifelse(TOLERANCE=='S',1,NA))
  inTaxa.1$MTOL <- with(inTaxa.1, ifelse(TOLERANCE=='I',1,NA))
  inTaxa.1$TOLR <- with(inTaxa.1, ifelse(TOLERANCE=='T',1,NA))

  # inTaxa.1 <- plyr::mutate(inTaxa,NTOL=ifelse(TOLERANCE %in% c('S','I'),1,NA)
  #                    ,INTL=ifelse(TOLERANCE=='S',1,NA)
  #                    ,MTOL=ifelse(TOLERANCE=='I',1,NA)
  #                    ,TOLR=ifelse(TOLERANCE=='T',1,NA)
  #                    )
  # Create empty data frames with all metric names in it
  empty_base <- data.frame(t(rep(NA,13)),stringsAsFactors=F)
  names(empty_base) <- c('INTLNTAX','INTLPIND','INTLPTAX',
                         'MTOLNTAX','MTOLPIND','MTOLPTAX','NTOLNTAX','NTOLPIND',
                         'NTOLPTAX','TOLRNTAX','TOLRPIND','TOLRPTAX','WTD_TV')

  # Add native metrics if CALCNAT='Y'
  if(CALCNAT=='Y'){
    empty_base.nat <- data.frame(t(rep(NA,17)),stringsAsFactors=F)
    names(empty_base.nat) <- c('NAT_INTLNTAX','NAT_INTLPIND','NAT_INTLPTAX',
                               'NAT_MTOLNTAX','NAT_MTOLPIND','NAT_MTOLPTAX','NAT_NTOLNTAX',
                               'NAT_NTOLPIND','NAT_NTOLPTAX','NAT_TOLRNTAX','NAT_TOLRPIND',
                               'NAT_TOLRPTAX','NAT_WTD_TV','NAT_PIND','NAT_PTAX','NAT_TOTLNIND','NAT_TOTLNTAX')
    empty_base <- cbind(empty_base,empty_base.nat)
  }


  if('VELOCITY' %in% names(inTaxa.1)){
    inTaxa.1$INTLRHEO <- with(inTaxa.1, ifelse(INTL==1 & VELOCITY=='R',1,NA))
    inTaxa.1$INTLLOT <- with(inTaxa.1, ifelse(INTL==1 & VELOCITY %in% c('R','O'),1,NA))
    # inTaxa.1 <- plyr::mutate(inTaxa.1, INTLRHEO=ifelse(INTL==1 & VELOCITY=='R',1,NA)
    #                    ,INTLLOT=ifelse(INTL==1 & VELOCITY %in% c('R','O'),1,NA)
    #                    )

    empty_vel <- data.frame(t(rep(NA,6)),stringsAsFactors=F)
    names(empty_vel) <- c('INTLRHEONTAX','INTLRHEOPIND','INTLRHEOPTAX',
                        'INTLLOTNTAX','INTLLOTPIND','INTLLOTPTAX')
    empty_base <- cbind(empty_base,empty_vel)

    if(CALCNAT=='Y'){
      empty_vel.nat <- empty_vel
      names(empty_vel.nat) <- paste('NAT',names(empty_vel),sep='_')
      empty_base <- cbind(empty_base,empty_vel.nat)
    }
  }

  if('HABITAT' %in% names(inTaxa.1)){
    inTaxa.1$NTOLBENT <- with(inTaxa.1, ifelse(NTOL==1 & HABITAT=='B',1,NA))
    # inTaxa.1 <- plyr::mutate(inTaxa.1, NTOLBENT=ifelse(NTOL==1 & HABITAT=='B',1,NA))

    empty_hab <- data.frame(t(rep(NA,3)),stringsAsFactors=F)
    names(empty_hab) <- c('NTOLBENTNTAX','NTOLBENTPIND','NTOLBENTPTAX')
    empty_base <- cbind(empty_base,empty_hab)

    if(CALCNAT=='Y'){
      empty_hab.nat <- empty_hab
      names(empty_hab.nat) <- paste('NAT',names(empty_hab),sep='_')
      empty_base <- cbind(empty_base,empty_hab.nat)
    }

  }

  if('TROPHIC' %in% names(inTaxa.1)){
    inTaxa.1$NTOLCARN <- with(inTaxa.1, ifelse(NTOL==1 & TROPHIC=='C',1,NA))
    inTaxa.1$INTLCARN <- with(inTaxa.1, ifelse(INTL==1 & TROPHIC=='C',1,NA))
    inTaxa.1$INTLINV <- with(inTaxa.1, ifelse(INTL==1 & TROPHIC=='I',1,NA))
    inTaxa.1$NTOLINV <- with(inTaxa.1, ifelse(NTOL==1 & TROPHIC=='I',1,NA))

    # inTaxa.1 <- plyr::mutate(inTaxa.1, NTOLCARN=ifelse(NTOL==1 & TROPHIC=='C',1,NA)
    #                    ,INTLCARN=ifelse(INTL==1 & TROPHIC=='C',1,NA)
    #                    ,INTLINV=ifelse(INTL==1 & TROPHIC=='I',1,NA)
    #                    ,NTOLINV=ifelse(NTOL==1 & TROPHIC=='I',1,NA))

    empty_trop <- data.frame(t(rep(NA,12)),stringsAsFactors=F)
    names(empty_trop) <- c('INTLINVNTAX','INTLINVPIND','INTLINVPTAX'
                           ,'INTLCARNNTAX','INTLCARNPIND','INTLCARNPTAX'
                           ,'NTOLCARNNTAX','NTOLCARNPIND','NTOLCARNPTAX'
                           ,'NTOLINVNTAX','NTOLINVPIND','NTOLINVPTAX')
    empty_base <- cbind(empty_base,empty_trop)

    if(CALCNAT=='Y'){
      empty_trop.nat <- empty_trop
      names(empty_trop.nat) <- paste('NAT',names(empty_trop),sep='_')
      empty_base <- cbind(empty_base,empty_trop.nat)
    }
  }

  if('MIGRATORY' %in% names(inTaxa.1)){
    inTaxa.1$INTLMIGR <- with(inTaxa.1, ifelse(INTL==1 & MIGRATORY=='Y',1,NA))
    # inTaxa.1 <- plyr::mutate(inTaxa.1,INTLMIGR=ifelse(INTL==1 & MIGRATORY=='Y',1,NA))

    empty_migr <- data.frame(t(rep(NA,3)),stringsAsFactors=F)
    names(empty_migr) <- c('INTLMIGRNTAX','INTLMIGRPIND','INTLMIGRPTAX')
    empty_base <- cbind(empty_base,empty_migr)

    if(CALCNAT=='Y'){
      empty_migr.nat <- empty_migr
      names(empty_migr.nat) <- paste('NAT',names(empty_migr),sep='_')
      empty_base <- cbind(empty_base,empty_migr.nat)
    }
  }

  params<-c('INTL','NTOL','MTOL','TOLR','INTLRHEO','INTLLOT','NTOLBENT','NTOLCARN','INTLCARN'
            ,'INTLINV','NTOLINV','INTLMIGR')

  inTaxa.2 <- subset(inTaxa.1,select=names(inTaxa.1) %in% c('TAXA_ID',params))

  params.use <- names(inTaxa.2)[names(inTaxa.2) %in% params]

  taxalong <- reshape(inTaxa.2, idvar = c('TAXA_ID'), direction = 'long',
                      varying = params.use, times = params.use, v.names = 'value', timevar = 'TRAIT')
  taxalong <- subset(taxalong, !is.na(value))

  # taxalong <- data.table::melt(inTaxa.2,id.vars='TAXA_ID',variable.name='TRAIT',na.rm=TRUE) %>%
  #   plyr::mutate(TRAIT=as.character(TRAIT))

  totals <- aggregate(x = list(TOTLNIND = inCts.1$TOTAL, TOTLNTAX = inCts.1$IS_DISTINCT), by = inCts.1[c('SAMPID')],
                      FUN = sum)

  inCts.2 <- merge(inCts.1, totals, by = 'SAMPID')
  inCts.2$CALCPIND <- with(inCts.2, TOTAL/TOTLNIND)
  inCts.2$CALCPTAX <- with(inCts.2, IS_DISTINCT/TOTLNTAX)

  # inCts.2 <- plyr::ddply(inCts.1, "SAMPID", mutate, TOTLNIND=sum(TOTAL),
  #                        TOTLNTAX=sum(IS_DISTINCT))

  if(CALCNAT=='Y'){
    inCts.2 <- inCts.2[,c('SAMPID','TOTAL','IS_DISTINCT','TAXA_ID','TOTLNTAX','TOTLNIND','NONNATIVE','CALCPIND','CALCPTAX')]
    # inCts.2 <- dplyr::select(inCts.2,SAMPID,TOTAL,IS_DISTINCT,TAXA_ID,TOTLNTAX,TOTLNIND,NONNATIVE)
  }else{
    inCts.2 <- inCts.2[,c('SAMPID','TOTAL','IS_DISTINCT','TAXA_ID','TOTLNTAX','TOTLNIND','CALCPIND','CALCPTAX')]
    # inCts.2 <- dplyr::select(inCts.2, SAMPID,TOTAL,IS_DISTINCT,TAXA_ID,TOTLNTAX,TOTLNIND)
  }

  # totals <- unique(inCts.2[,c('SAMPID','TOTLNTAX','TOTLNIND')])

  # Merge the count data with the taxalist containing only the traits of
  # interest
  traitDF <- merge(inCts.2, taxalong, by='TAXA_ID')

  # Calculate no. individuals, % individuals, no. taxa, and % taxa for each
  # trait in taxalist
  outMet.1 <- aggregate(x = list(NTAX = traitDF$IS_DISTINCT),
                        by = traitDF[c('SAMPID','TRAIT')],
                        FUN = sum)

  outMet.2 <- aggregate(x = list(PIND = traitDF$CALCPIND, PTAX = traitDF$CALCPTAX),
                        by = traitDF[c('SAMPID','TRAIT')],
                        FUN = function(x){round(sum(x)*100, 2)})

  outMet <- merge(outMet.1, outMet.2, by = c('SAMPID','TRAIT'))


  # outMet <- plyr::ddply(traitDF, c("SAMPID", "TRAIT"), summarise,
  #                 NTAX=sum(IS_DISTINCT),
  #                 PIND=round(sum(TOTAL/TOTLNIND)*100,2),
  #                 PTAX=round(sum(IS_DISTINCT/TOTLNTAX)*100,2), .progress='tk')

  # Melt df to create metric names, then recast into wide format with metric
  # names
  outLong <- reshape(outMet, idvar = c('SAMPID','TRAIT'), direction = 'long',
                     varying = c('NTAX','PIND','PTAX'), timevar = 'variable',
                     v.names = 'value', times = c('NTAX','PIND','PTAX'))

  # outLong <- data.table::melt(outMet,id.vars=c('SAMPID','TRAIT'))
  outLong$variable <- paste(outLong$TRAIT,outLong$variable,sep='')
  outLong$TRAIT <- NULL

  outWide <- reshape(outLong, idvar = 'SAMPID', direction = 'wide',
                     v.names = 'value', timevar = 'variable')
  names(outWide) <- gsub("value\\.", "", names(outWide))

  outWide <- merge(outWide, totals, by = 'SAMPID', all.y=TRUE)

  # outWide <- data.table::dcast(outLong,SAMPID~variable,value.var='value') %>%
  #   merge(totals,by='SAMPID',all.y=T)

  if(nrow(subset(inTaxa,!is.na(TOL_VAL)))>0){
    TVI <- tolindexFish(inCts.2,inTaxa)
    outWide <- merge(outWide,TVI,by="SAMPID",all.x=TRUE)
  }

  # Now run native metrics if CALCNAT='Y'
  ## If the variable NON_NATIVE is included and populated in inCts1, create inNative data frame
  if(CALCNAT=='Y'){
    if(any(unique(inCts.2$NON_NATIVE) %nin% c('Y','N'))){
      return(print("No native and alien datasets were created because NON_NATIVE must only be 'Y' or 'N' values"))
    }else{
      inNative <- subset(inCts.2,NONNATIVE=='N')
      if(length(inNative)>0){
        natTot <- aggregate(x = list(NAT_TOTLNIND = inNative$TOTAL, NAT_TOTLNTAX = inNative$IS_DISTINCT),
                            by = inNative[c('SAMPID')], FUN = sum)

        inNative.tot <- merge(inNative, natTot, by = 'SAMPID')
        # inNative.tot <- plyr::ddply(inNative,c('SAMPID'),mutate,NAT_TOTLNIND=sum(TOTAL),
        #                                      NAT_TOTLNTAX=sum(IS_DISTINCT))
        # totals.nat <- unique(inNative.tot[,c('SAMPID','NAT_TOTLNIND','NAT_TOTLNTAX')])

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

        # natMets <- merge(inNative.tot, taxalong, by='TAXA_ID') %>%
        #   plyr::ddply(c('SAMPID','TRAIT','NAT_TOTLNTAX','NAT_TOTLNIND'),summarise,
        #                           NTAX=sum(IS_DISTINCT),
        #                           PIND=round(sum(TOTAL/NAT_TOTLNIND)*100,2),
        #                           PTAX=round(sum(IS_DISTINCT/NAT_TOTLNTAX)*100,2), .progress='tk')

        natMets.long <- reshape(natMets.comb, idvar = c('SAMPID','TRAIT'), direction = 'long',
                                varying = c('NTAX', 'PIND', 'PTAX'), timevar = 'variable', v.names = 'value',
                                times = c('NTAX', 'PIND', 'PTAX'))
        natMets.long$variable <- with(natMets.long, paste('NAT_',TRAIT,variable,sep=''))
        natMets.long$TRAIT <- NULL

        # natMets.long <- data.table::melt(natMets,id.vars=c('SAMPID','TRAIT','NAT_TOTLNTAX','NAT_TOTLNIND')) %>%
        #   mutate(variable=paste('NAT_',TRAIT,variable,sep=''))

        natMets.fin <- reshape(natMets.long, idvar = c('SAMPID'), direction = 'wide',
                               v.names = 'value', timevar = 'variable')
        names(natMets.fin) <- gsub("value\\.", "", names(natMets.fin))
        natMets.fin <- merge(natMets.fin, natTot, by = 'SAMPID', all.y=TRUE)

        # natMets.1 <- data.table::dcast(natMets.long,SAMPID~variable,value.var='value') %>%
        #   merge(totals.nat,by='SAMPID',all.y=T)

        outWide.1 <- merge(outWide, natMets.fin, all = TRUE)
        # outWide.1 <- merge(outWide,natMets.1,all=T)

        if(nrow(subset(inTaxa,!is.na(TOL_VAL)))>0){
          TVI <- tolindexFish(inNative,inTaxa)
          names(TVI)[names(TVI)=='WTD_TV'] <- 'NAT_WTD_TV'
          # TVI <- tolindexFish(inNative,inTaxa) %>%
          #   plyr::rename(c('WTD_TV'='NAT_WTD_TV'))
          outWide.1 <- merge(outWide.1,TVI,by="SAMPID",all.x=TRUE)
        }

        outWide.1$NAT_PTAX <- with(outWide.1, round((NAT_TOTLNTAX/TOTLNTAX)*100,2))
        outWide.1$NAT_PIND <- with(outWide.1, round((NAT_TOTLNIND/TOTLNIND)*100,2))
        outWide.1$TOTLNTAX <- NULL
        outWide.1$TOTLNIND <- NULL
        # outWide.1 <- plyr::mutate(outWide.1,NAT_PTAX=round((NAT_TOTLNTAX/TOTLNTAX)*100,2),NAT_PIND=round((NAT_TOTLNIND/TOTLNIND)*100,2)) %>%
        #   select(-TOTLNTAX,-TOTLNIND)

      }else{
        outWide.1 <- outWide
        outWide.1$TOTLNTAX <- NULL
        outWide.1$TOTLNIND <- NULL
        # outWide.1 <- select(outWide,-TOTLNTAX,-TOTLNIND)
      }
    }
    }else{
      outWide.1 <- outWide
      outWide.1$TOTLNTAX <- NULL
      outWide.1$TOTLNIND <- NULL
      # outWide.1 <- select(outWide,-TOTLNTAX,-TOTLNIND)
      }

  outWide.all <- merge(outWide.1, subset(empty_base), all = TRUE)
  outwide.all <- outWide.all[!is.na(outWide.all$SAMPID),]
  outWide.all <- merge(outWide.all, samples, by = 'SAMPID', all.y = TRUE)

  # outWide.all <- merge(outWide.1,empty_base,all=TRUE) %>%
  #   filter(!is.na(SAMPID)) %>%
  #   merge(samples,by='SAMPID',all.y=T)

  # If we re-melt df now, we have missing values where the metric should be a
  # zero, so we can set NAs to 0 now
  # outWide.all[is.na(outWide.all)] <- 0
  updNames <- names(outWide.all)[names(outWide.all) %nin% c('WTD_TV','NAT_WTD_TV','SAMPID')]
  outWide.all[,updNames] <- lapply(outWide.all[,updNames], function(x){ifelse(is.na(x), 0, x)})
  # outLong.1 <- data.table::melt(outWide.all,id.vars=c(sampID,'SAMPID')) %>%
  #   plyr::mutate(value=ifelse(is.na(value) & variable %nin% c('WTD_TV','NAT_WTD_TV'),0,value))
  #
  # # Finally, we can recast the metrics df into wide format for output
  # lside <- paste(paste(sampID,collapse='+'),'SAMPID',sep='+')
  # formula <- paste(lside,'~variable',sep='')
  # outWide.2 <- data.table::dcast(outLong.1,eval(formula),value.var='value')

  # Merge metrics with the original indata so that those without metrics because
  # no sample was collected are still output with missing values
  outAll <- merge(outWide.all, totals, by = 'SAMPID', all.x = TRUE)
  outAll$SAMPID <- NULL
  # outAll <- merge(outWide.2,totals,by='SAMPID',all.x=T) %>%
  #   dplyr::select(-SAMPID)

  return(outAll)

}
