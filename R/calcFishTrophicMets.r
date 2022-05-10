#' @export
#' @title Calculate all trophic-related metrics
#' @description This function calculates all of the trophic metrics,
#' and if habitat values are included in the taxalist, additional
#' metrics are calculated.
#' @param indata Input data frame containing variables as identified
#' in the arguments for \emph{sampID}, \emph{dist}, \emph{ct},
#' \emph{taxa_id}, as well as the optional variable
#' for non-native status in \emph{nonnat}.
#' @param inTaxa Data frame containing fish taxalist, along with autecology
#' traits. At a minimum, this taxalist must contain variables matching
#' the argument for \emph{trophic}. If habitat values are included, as
#' specified by the argument \emph{habitat} additional metrics that combine
#' trophic status and habitat are also calculated. The variable specified
#' in the argument \emph{taxa_id} must match that in \emph{indata}.
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
#' @param habitat An optional string with the name of the variable in
#' \emph{inTaxa} containing habitat preference values. Valid
#' values include B (benthic), W (water column), E (edge), or
#' blank, if unknown. The default value is \emph{HABITAT}.
#' @param trophic A string with the name of the variable in
#' \emph{inTaxa} containing trophic category values. Valid
#' values include I (invertivore), C (carnivore), O (omnivore),
#' H (herbivore), or blank if unknown. The default value is
#' \emph{inTaxa} is TROPHIC.
#' @param nonnat A string with the name of the optional variable in
#' \emph{inCts} containing non-native status. Valid values are 'Y' for
#' non-native and 'N' for native. The default name
#' is \emph{NONNATIVE}.
#' @return A data frame containing the variables in sampID and
#' the fish trophic metrics as additional variables. Metric
#' descriptions are included in \emph{NRSA_Fish_Metric_Descriptions.pdf},
#' included in this package. The names of
#' metrics include  INVNTAX, INVPIND, INVPTAX, CARNNTAX, CARNPIND,
#' CARNPTAX, OMNINTAX, OMNIPIND, OMNIPTAX, HERBNTAX, HERBPIND,
#' HERBPTAX, TOTLNIND, and TOTLNTAX.
#'
#' If a non-native status variable is included, these metrics are also
#' calculated:
#' NAT_INVNTAX, NAT_INVPIND, NAT_INVPTAX, NAT_CARNNTAX,
#' NAT_CARNPIND, NAT_CARNPTAX, NAT_OMNINTAX, NAT_OMNIPIND,
#' NAT_OMNIPTAX, NAT_HERBNTAX, NAT_HERBPIND, NAT_HERBPTAX,
#' NAT_TOTLNTAX, NAT_TOTLNIND, NAT_PIND, NAT_PTAX.
#'
#' Additional metrics calculated if habitat is included in
#' \emph{inTaxa}: BENTINVNIND, BENTINVNTAX, BENTINVPIND.
#'
#' If a non-native status variable is included, these additional metrics
#' are also calculated: NAT_BENTINVNIND, NAT_BENTINVNTAX, NAT_BENTINVPIND.
#'
#' @author Karen Blocksom \email{Blocksom.Karen@epa.gov}
#'
calcFishTrophicMets <- function(indata, inTaxa=NULL, sampID='UID', dist='IS_DISTINCT'
                            , ct='TOTAL', taxa_id='TAXA_ID', habitat='HABITAT'
                            , trophic='TROPHIC', nonnat='NONNATIVE'){

  # Convert data into data frames just in case
  indata <- as.data.frame(indata)
  inTaxa <- as.data.frame(inTaxa)

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
  necTraits <- c(trophic)
  if(any(necTraits %nin% names(inTaxa))){
    msgTraits <- which(necTraits %nin% names(inTaxa))
    return(paste("Some of the traits are missing from the taxa list. The following are \nrequired for metric calculations to run:", necTraits[msgTraits]))
  }
  optTraits <- c(habitat)
  if(any(optTraits %nin% names(inTaxa))){
    msgTraits <- which(optTraits %nin% names(inTaxa))
    print(paste("Optional trait is missing from the taxa list. Any tolerance metrics also using this trait will not be calculated:\n",
                optTraits[msgTraits],"\n"))
  }


  inTaxa <- subset(inTaxa,select=names(inTaxa) %in% c(taxa_id,habitat,trophic))

  # Rename counts and distinct variables to TOTAL and IS_DISTINCT
  names(indata)[names(indata)==ct] <- 'TOTAL'
  names(indata)[names(indata)==dist] <- 'IS_DISTINCT'
  names(indata)[names(indata)==taxa_id] <- 'TAXA_ID'
  names(indata)[names(indata)==nonnat] <- 'NONNATIVE'
  names(inTaxa)[names(inTaxa)==taxa_id] <- 'TAXA_ID'
  names(inTaxa)[names(inTaxa)==trophic] <- 'TROPHIC'

  if(habitat %in% names(inTaxa)){
    names(inTaxa)[names(inTaxa)==habitat] <- 'HABITAT'
  }

  indata[,c('TOTAL','IS_DISTINCT')] <- lapply(indata[,c('TOTAL','IS_DISTINCT')],as.numeric)
  indata$TAXA_ID <- as.character(indata$TAXA_ID)
  inTaxa$TAXA_ID <- as.character(inTaxa$TAXA_ID)

  ## for inCts1, keep only observations without missing or zero TOTAL values or TAXA_ID and TAXA_ID!=99999
  indata.1 <- subset(indata,!is.na(TAXA_ID) & !is.na(TOTAL) & TOTAL!=0)

  ## Now sum by TAXA_ID for ANOM_CT and for TOTAL for each sample
  # Two approaches depending on whether or not NON_NATIVE occurs in the counts data frame
  if('NONNATIVE' %in% names(indata.1)){
    maxDist <- aggregate(x = list(IS_DISTINCT = indata.1$IS_DISTINCT), by = indata.1[c('SAMPID','TAXA_ID','NONNATIVE')],
                         FUN = function(x){max(as.integer(x))})
    sumTot <- aggregate(x = list(TOTAL = indata.1$TOTAL), by = indata.1[c('SAMPID','TAXA_ID','NONNATIVE')], FUN = sum)

    indata.2 <- merge(maxDist, sumTot, by = c('SAMPID','TAXA_ID','NONNATIVE'))
    CALCNAT <- 'Y'
  }else{
    maxDist <- aggregate(x = list(IS_DISTINCT = indata.1$IS_DISTINCT), by = indata.1[c('SAMPID','TAXA_ID')],
                         FUN = function(x){max(as.integer(x))})
    sumTot <- aggregate(x = list(TOTAL = indata.1$TOTAL), by = indata.1[c('SAMPID','TAXA_ID')], FUN = sum)

    indata.2 <- merge(maxDist, sumTot, by = c('SAMPID','TAXA_ID'))
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

  inCts.1 <- inCts[inCts$TAXA_ID %in% inTaxa$TAXA_ID,]

  if(CALCNAT=='Y'){
    inCts.1 <- inCts.1[!is.na(inCts.1$TOTAL) & inCts.1$TOTAL>0,c('SAMPID','TAXA_ID','TOTAL','IS_DISTINCT','NONNATIVE')]
  }else{
    inCts.1 <- inCts.1[!is.na(inCts.1$TOTAL) & inCts.1$TOTAL>0,c('SAMPID','TAXA_ID','TOTAL','IS_DISTINCT')]
  }
  # Now create indicator variables
  inTaxa.1 <- inTaxa
  inTaxa.1$INV <- with(inTaxa.1, ifelse(TROPHIC %in% c('I'),1,NA))
  inTaxa.1$CARN <- with(inTaxa.1, ifelse(TROPHIC=='C',1,NA))
  inTaxa.1$OMNI <- with(inTaxa.1, ifelse(TROPHIC=='O',1,NA))
  inTaxa.1$HERB <- with(inTaxa.1, ifelse(TROPHIC=='H',1,NA))

  # Create empty data frames with all metric names in it
  empty_base <- data.frame(t(rep(NA,12)),stringsAsFactors=F)
  names(empty_base) <- c('INVNTAX','INVPIND','INVPTAX',
                         'CARNNTAX','CARNPIND','CARNPTAX',
                         'OMNINTAX','OMNIPIND','OMNIPTAX',
                         'HERBNTAX','HERBPIND','HERBPTAX')

  # Add native metrics if CALCNAT='Y'
  if(CALCNAT=='Y'){
    empty_base.nat <- empty_base
    names(empty_base.nat) <- paste('NAT',names(empty_base),sep='_')
    empty_base <- cbind(empty_base,empty_base.nat)
  }

  if('HABITAT' %in% names(inTaxa.1)){
    inTaxa.1$BENTINV <- with(inTaxa.1, ifelse(HABITAT=='B' & TROPHIC=='I',1,NA))

    empty_hab <- data.frame(t(rep(NA,3)),stringsAsFactors=F)
    names(empty_hab) <- c('BENTINVNTAX','BENTINVPIND','BENTINVPTAX')
    empty_base <- cbind(empty_base,empty_hab)

    if(CALCNAT=='Y'){
      empty_hab.nat <- empty_hab
      names(empty_hab.nat) <- paste('NAT',names(empty_hab),sep='_')
      empty_base <- cbind(empty_base,empty_hab.nat)
    }

  }

  params<-c('INV','CARN','OMNI','HERB','BENTINV')

  inTaxa.2 <- subset(inTaxa.1,select=names(inTaxa.1) %in% c('TAXA_ID',params))

  params.use <- names(inTaxa.2)[names(inTaxa.2) %in% params]

  taxalong <- reshape(inTaxa.2, idvar = c('TAXA_ID'), direction = 'long',
                      varying = params.use, times = params.use, v.names = 'value', timevar = 'TRAIT')
  taxalong <- subset(taxalong, !is.na(value))

  totals <- aggregate(x = list(TOTLNIND = inCts.1$TOTAL, TOTLNTAX = inCts.1$IS_DISTINCT), by = inCts.1[c('SAMPID')],
                      FUN = sum)

  inCts.2 <- merge(inCts.1, totals, by = 'SAMPID')
  inCts.2$CALCPIND <- with(inCts.2, TOTAL/TOTLNIND)
  inCts.2$CALCPTAX <- with(inCts.2, IS_DISTINCT/TOTLNTAX)

  if(CALCNAT=='Y'){
    inCts.2 <- inCts.2[,c('SAMPID','TOTAL','IS_DISTINCT','TAXA_ID','TOTLNTAX','TOTLNIND','NONNATIVE','CALCPIND','CALCPTAX')]
  }else{
    inCts.2 <- inCts.2[,c('SAMPID','TOTAL','IS_DISTINCT','TAXA_ID','TOTLNTAX','TOTLNIND','CALCPIND','CALCPTAX')]
  }

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

        natMets.fin <- merge(natMets.fin, natTot, by = 'SAMPID', all.y=TRUE)

        outWide.1 <- merge(outWide, natMets.fin, all = TRUE)

        outWide.1$NAT_PTAX <- with(outWide.1, round((NAT_TOTLNTAX/TOTLNTAX)*100,2))
        outWide.1$NAT_PIND <- with(outWide.1, round((NAT_TOTLNIND/TOTLNIND)*100,2))
        outWide.1$TOTLNTAX <- NULL
        outWide.1$TOTLNIND <- NULL

      }else{
        outWide.1 <- outWide
        outWide.1$TOTLNTAX <- NULL
        outWide.1$TOTLNIND <- NULL
      }
    }
  }else{
    outWide.1 <- outWide
    outWide.1$TOTLNTAX <- NULL
    outWide.1$TOTLNIND <- NULL
  }

  outWide.all <- merge(outWide.1, subset(empty_base), all = TRUE)
  outwide.all <- outWide.all[!is.na(outWide.all$SAMPID),]
  outWide.all <- merge(outWide.all, samples, by = 'SAMPID', all.y = TRUE)

  # If we re-melt df now, we have missing values where the metric should be a
  # zero, so we can set NAs to 0 now
  outWide.all[is.na(outWide.all)] <- 0

  # Merge metrics with the original indata so that those without metrics because
  # no sample was collected are still output with missing values
  outAll <- merge(outWide.all, totals, by = 'SAMPID', all.x = TRUE)
  outAll$SAMPID <- NULL

  return(outAll)

}
