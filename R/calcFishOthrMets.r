#' @export
#' @title Calculate other subsets of fish metrics
#' @description This function calculates other metrics using the
#' traits for velocity preference, migratory tendency, reproductive habits,
#' and temperature preference, depending on what traits are included
#' in the input taxalist. Native status versions of metrics are calculated
#' if a non-native status variable is included in the input count data.
#' @param indata Input data frame containing variables as identified
#' in the arguments for \emph{sampID}, \emph{dist}, \emph{ct},
#' \emph{taxa_id}, as well as the optional variable
#' for non-native status in \emph{nonnat}.
#' @param inTaxa Data frame containing fish taxalist, along with autecology
#' traits. At a minimum, this taxalist must contain at least one of the
#' variables matching \emph{migr}, \emph{vel}, \emph{reprod}, \emph{temp}.
#' The variable specified in the argument \emph{taxa_id}
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
#' @param vel A strings with the name of the variable in the
#' \emph{inTaxa} taxalist containing velocity preference values.
#' Valid values include R (Rheophil), P (Pool), O (Other), or
#' blank, if unknown. The default value is \emph{VELOCITY}.
#' @param migr A string with the name of the variable in
#' \emph{inTaxa} containing migratory status value. Valid
#' values include N (No), Y (Yes), and blank if unknown. The
#' default value is \emph{MIGRATORY}.
#' @param reprod A string with the name of the variable in
#' \emph{inTaxa} containing reproductive trait values. Valid
#' values include C (Clean, coarse lithophil), D (Drifter),
#' G (Guarder), O (Other), or blank if unknown. The default
#' value is \emph{REPROD}.
#' @param temp A string with the name of the variable in
#' \emph{inTaxa} containing the temperature preference values.
#' Valid values include WM (warm), CD (cold water), CL (cool water),
#' or blank if unknow. The default value is \emph{TEMP}.
#' @param nonnat A string with the name of the optional variable in
#' \emph{inCts} containing non-native status. Valid values are 'Y' for
#' non-native and 'N' for native. The default name
#' is \emph{NONNATIVE}.
#' @return A data frame containing the variables in sampID and
#' the fish metrics as additional variables. Metric
#' descriptions are included in \emph{NRSA_Fish_Metric_Descriptions.pdf},
#' included in this package. The names of
#' metrics include  RHEONTAX, RHEOPIND, RHEOPTAX, LOTNTAX, LOTPIND,
#' LOTPTAX,  MIGRNTAX, MIGRPIND, MIGRPTAX,  LITHNTAX, LITHPIND,
#' LITHPTAX,  COLDNTAX, COLDPIND, COLDPTAX, TOTLNIND,
#' and TOTLNTAX.
#'
#' If a non-native status variable is included, these metrics are also
#' calculated:
#' NAT_RHEONTAX, NAT_RHEOPIND, NAT_RHEOPTAX, NAT_LOTNTAX,
#' NAT_LOTPIND, NAT_LOTPTAX, NAT_ MIGRNTAX, NAT_MIGRPIND, NAT_MIGRPTAX,
#' NAT_ LITHNTAX, NAT_LITHPIND, NAT_LITHPTAX, NAT_ COLDNTAX, NAT_COLDPIND,
#' NAT_COLDPTAX, NAT_NAT_TOTLNTAX, NAT_NAT_TOTLNIND, NAT_PIND, NAT_PTAX.
#'
#'
#' @author Karen Blocksom \email{Blocksom.Karen@epa.gov}
#'
calcFishOtherMets <- function(indata, inTaxa=NULL, sampID='UID', dist='IS_DISTINCT'
                            , ct='TOTAL', taxa_id='TAXA_ID', vel='VELOCITY'
                            , migr='MIGRATORY', reprod='REPROD', temp='TEMP'
                            , nonnat='NONNATIVE'){
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
  optTraits <- c(migr,vel,reprod,temp)
  if(all(optTraits %nin% names(inTaxa))){
    return(print("Need at least one of the following in inTaxa taxalist: migr, vel, reprod, temp"))
  }

  if(any(optTraits %nin% names(inTaxa))){
    msgTraits <- which(optTraits %nin% names(inTaxa))
    print(paste("Some optional traits are missing from the taxa list. Any tolerance metrics also using these traits will not be calculated:\n",
                optTraits[msgTraits],"\n"))
  }


  inTaxa <- subset(inTaxa,select=names(inTaxa) %in% c(taxa_id,vel,migr,reprod,temp))

  # Rename counts and distinct variables to TOTAL and IS_DISTINCT
  names(indata)[names(indata)==ct] <- 'TOTAL'
  names(indata)[names(indata)==dist] <- 'IS_DISTINCT'
  names(indata)[names(indata)==taxa_id] <- 'TAXA_ID'
  names(indata)[names(indata)==nonnat] <- 'NONNATIVE'
  names(inTaxa)[names(inTaxa)==taxa_id] <- 'TAXA_ID'

  if(vel %in% names(inTaxa)){
    names(inTaxa)[names(inTaxa)==vel] <- 'VELOCITY'
  }

  if(migr %in% names(inTaxa)){
    names(inTaxa)[names(inTaxa)==migr] <- 'MIGRATORY'
  }

  if(reprod %in% names(inTaxa)){
    names(inTaxa)[names(inTaxa)==reprod] <- 'REPROD'
  }

  if(temp %in% names(inTaxa)){
    names(inTaxa)[names(inTaxa)==temp] <- 'TEMP'
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

  # Create empty data frames with all metric names in it
  empty_base <- data.frame(SAMPID=NA,stringsAsFactors=F)

  if('VELOCITY' %in% names(inTaxa.1)){
    inTaxa.1$RHEO <- with(inTaxa.1, ifelse(VELOCITY=='R',1,NA))
    inTaxa.1$LOT <- with(inTaxa.1, ifelse(VELOCITY %in% c('R','O'),1,NA))

    empty_vel <- data.frame(t(rep(NA,6)),stringsAsFactors=F)
    names(empty_vel) <- c('RHEONTAX','RHEOPIND','RHEOPTAX',
                          'LOTNTAX','LOTPIND','LOTPTAX')
    empty_base <- cbind(empty_base,empty_vel)

    if(CALCNAT=='Y'){
      empty_vel.nat <- empty_vel
      names(empty_vel.nat) <- paste('NAT',names(empty_vel),sep='_')
      empty_base <- cbind(empty_base,empty_vel.nat)
    }
  }

  if('MIGRATORY' %in% names(inTaxa.1)){
    inTaxa.1$MIGR <- with(inTaxa.1, ifelse(MIGRATORY=='Y',1,NA))

    empty_migr <- data.frame(t(rep(NA,3)),stringsAsFactors=F)
    names(empty_migr) <- c('MIGRNTAX','MIGRPIND','MIGRPTAX')
    empty_base <- cbind(empty_base,empty_migr)

    if(CALCNAT=='Y'){
      empty_migr.nat <- empty_migr
      names(empty_migr.nat) <- paste('NAT',names(empty_migr),sep='_')
      empty_base <- cbind(empty_base,empty_migr.nat)
    }
  }

  if('REPROD' %in% names(inTaxa.1)){
    inTaxa.1$LITH <- with(inTaxa.1, ifelse(REPROD=='C',1,NA))

    empty_repr <- data.frame(t(rep(NA,3)),stringsAsFactors=F)
    names(empty_repr) <- c('LITHNTAX','LITHPIND','LITHPTAX')
    empty_base <- cbind(empty_base,empty_repr)

    if(CALCNAT=='Y'){
      empty_repr.nat <- empty_repr
      names(empty_repr.nat) <- paste('NAT',names(empty_repr),sep='_')
      empty_base <- cbind(empty_base,empty_repr.nat)
    }
  }

  if('TEMP' %in% names(inTaxa.1)){
    inTaxa.1$COLD <- with(inTaxa.1, ifelse(TEMP=='CD',1,NA))

    empty_temp <- data.frame(t(rep(NA,3)),stringsAsFactors=F)
    names(empty_temp) <- c('COLDNTAX','COLDPIND','COLDPTAX')
    empty_base <- cbind(empty_base,empty_temp)

    if(CALCNAT=='Y'){
      empty_temp.nat <- empty_temp
      names(empty_temp.nat) <- paste('NAT',names(empty_temp),sep='_')
      empty_base <- cbind(empty_base,empty_temp.nat)
    }
  }

  params<-c('COLD','LITH','RHEO','LOT','MIGR')

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

  outWide.all <- merge(outWide.1, empty_base, all = TRUE)
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
