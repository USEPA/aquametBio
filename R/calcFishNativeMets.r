#' @export
#' @title Calculate all native status metrics
#' @description This function calculates all of the metrics based
#' only on native status, including alien metrics.
#' @param indata Input data frame containing variables as identified
#' in the arguments for \emph{sampID}, \emph{dist}, \emph{ct},
#' \emph{taxa_id}, and \emph{nonnat}.
#' @param sampID sampID A character vector containing the names of all
#' variables in indf that specify a unique sample. If not specified,
#' the default is \emph{UID}.
#' @param dist A string with the name of the distinctness variable,
#' which is assumed to have only values of 0 or 1. If not specified,
#' the default is \emph{IS_DISTINCT}.
#' @param ct A string with the name of the count variable. If not
#' specified, the default is \emph{TOTAL}.
#' @param taxa_id A string with the name of the taxon variable
#' in \emph{indata}. The default value is \emph{TAXA_ID}.
#' @param nonnat A string with the name of the optional variable in
#' \emph{inCts} containing non-native status. Valid values are 'Y' for
#' non-native and 'N' for native. The default name
#' is \emph{NONNATIVE}.
#' @return A data frame containing the variables in sampID and
#' the fish native status metrics as additional variables. Metric
#' descriptions are included in \emph{NRSA_Fish_Metric_Descriptions.pdf},
#' included in this package. The names of
#' metrics include  NAT_TOTLNTAX, NAT_PTAX, NAT_PIND, ALIENNTAX,
#' ALIENPTAX, ALIENPIND.
#'
#' @author Karen Blocksom \email{Blocksom.Karen@epa.gov}
#'
calcFishNativeMets <- function(indata, sampID='UID', dist='IS_DISTINCT'
                              , ct='TOTAL', taxa_id='TAXA_ID', nonnat='NONNATIVE'){

  ctVars <- c(sampID,dist,ct,taxa_id,nonnat)
  if(any(ctVars %nin% names(indata))){
    msgTraits <- which(ctVars %nin% names(indata))
    print(paste("Missing variables in input data frame:",paste(ctVars[msgTraits],collapse=',')))
    return(NULL)
  }

  # Combine all values in sampID into one sampID in df
  for(i in 1:length(sampID)){
    if(i==1) indata$SAMPID <- indata[,sampID[i]]
    else indata$SAMPID <- paste(indata$SAMPID,indata[,sampID[i]],sep='.')
  }
  # Keep data frame with crosswalk info between sampID and SAMPID
  samples <- unique(subset(indata,select=c(sampID,'SAMPID')))

  # Rename counts and distinct variables to TOTAL and IS_DISTINCT
  names(indata)[names(indata)==ct] <- 'TOTAL'
  names(indata)[names(indata)==dist] <- 'IS_DISTINCT'
  names(indata)[names(indata)==taxa_id] <- 'TAXA_ID'
  names(indata)[names(indata)==nonnat] <- 'NONNATIVE'

  indata[,c('TOTAL','IS_DISTINCT')] <- lapply(indata[,c('TOTAL','IS_DISTINCT')],as.numeric)
  indata$TAXA_ID <- as.character(indata$TAXA_ID)

  ## for inCts1, keep only observations without missing or zero TOTAL values or TAXA_ID and TAXA_ID!=99999
  indata.1 <- subset(indata,!is.na(TAXA_ID) & !is.na(TOTAL) & TOTAL!=0)

  maxDist <- aggregate(x = list(IS_DISTINCT = indata.1$IS_DISTINCT), by = indata.1[c('SAMPID','TAXA_ID','NONNATIVE')],
                       FUN = function(x){max(as.integer(x))})
  sumTot <- aggregate(x = list(TOTAL = indata.1$TOTAL), by = indata.1[c('SAMPID','TAXA_ID','NONNATIVE')], FUN = sum)

  inCts <- merge(maxDist, sumTot, by = c('SAMPID','TAXA_ID','NONNATIVE'))
#  indata.2 <- plyr::ddply(indata.1,c('SAMPID','TAXA_ID','NONNATIVE'),summarise,IS_DISTINCT=max(IS_DISTINCT),TOTAL=sum(TOTAL))

  # Make sure all necessary columns in inCts are numeric
#  inCts <- plyr::mutate(indata.2,TOTAL=as.numeric(TOTAL),IS_DISTINCT=as.integer(IS_DISTINCT))

  if(any(unique(inCts$NON_NATIVE) %nin% c('Y','N'))){
    return(print("No native and alien datasets were created because NON_NATIVE must only be 'Y' or 'N' values"))
  }

  totals <- aggregate(x = list(TOTLNIND = inCts$TOTAL, TOTLNTAX = inCts$IS_DISTINCT), by = inCts[c('SAMPID')],
                      FUN = sum)

  inCts.1 <- merge(inCts, totals, by = 'SAMPID')
  inCts.1$ALIEN <- with(inCts.1, ifelse(NONNATIVE=='Y',1,NA))
  inCts.1$NAT <- with(inCts.1, ifelse(NONNATIVE=='N',1,NA))
  # inCts.1 <- plyr::ddply(inCts, "SAMPID", mutate, TOTLNIND=sum(TOTAL),TOTLNTAX=sum(IS_DISTINCT)) %>%
  #   mutate(ALIEN=ifelse(NONNATIVE=='Y',1,NA),NAT=ifelse(NONNATIVE=='N',1,NA))

  # totals <- unique(inCts.1[,c('SAMPID','TOTLNTAX','TOTLNIND')])

  empty_base <- data.frame(t(rep(NA,7)),stringsAsFactors=F)
  names(empty_base) <- c('NAT_TOTLNIND','NAT_TOTLNTAX','NAT_PIND','NAT_PTAX','ALIENPIND','ALIENNTAX','ALIENPTAX')

  inCts.1$alienntax <- with(inCts.1, IS_DISTINCT*ALIEN)
  inCts.1$alienpind <- with(inCts.1, TOTAL*ALIEN/TOTLNIND)
  inCts.1$alienptax <- with(inCts.1, IS_DISTINCT*ALIEN/TOTLNTAX)
  inCts.1$natntax <- with(inCts.1, IS_DISTINCT*NAT)
  inCts.1$natnind <- with(inCts.1, TOTAL*NAT)
  inCts.1$natpind <- with(inCts.1, TOTAL*NAT/TOTLNIND)
  inCts.1$natptax <- with(inCts.1, IS_DISTINCT*NAT/TOTLNTAX)

  outMet.1 <- with(inCts.1, aggregate(x = list(ALIENNTAX = alienntax, NAT_TOTLNTAX = natntax, NAT_TOTLNIND = natnind),
                                      by = inCts.1[c('SAMPID')], FUN = function(x){sum(x, na.rm=T)}))
  outMet.2 <- with(inCts.1, aggregate(x = list(ALIENPIND = alienpind, ALIENPTAX = alienptax, NAT_PIND = natpind, NAT_PTAX = natptax),
                                      by = inCts.1[c('SAMPID')], FUN = function(x){round(sum(x, na.rm=T)*100, 2)}))

  outMet <- merge(outMet.1, outMet.2, by = 'SAMPID')
  # outMet <- plyr::ddply(inCts.1, c("SAMPID"), summarise,
  #                       ALIENNTAX=sum(IS_DISTINCT*ALIEN,na.rm=T),
  #                       ALIENPIND=round(sum(TOTAL*ALIEN/TOTLNIND,na.rm=T)*100,2),
  #                       ALIENPTAX=round(sum(IS_DISTINCT*ALIEN/TOTLNTAX,na.rm=T)*100,2),
  #                       NAT_NTAX=sum(IS_DISTINCT*NAT,na.rm=T),
  #                       NAT_NIND=sum(TOTAL*NAT,na.rm=T),
  #                       NAT_PIND=round(sum(TOTAL*NAT/TOTLNIND,na.rm=T)*100,2),
  #                       NAT_PTAX=round(sum(IS_DISTINCT*NAT/TOTLNTAX,na.rm=T)*100,2),.progress='tk')

  outMet.3 <- merge(outMet,empty_base,all=TRUE)
  outMet.3 <- outMet.3[!is.na(outMet.3$SAMPID),]
  outMet.3 <- merge(outMet.3, samples, by = 'SAMPID', all.y=T)

  # outMet.1 <- merge(outMet,empty_base,all=TRUE) %>%
  #   dplyr::filter(!is.na(SAMPID)) %>%
  #   merge(samples,by='SAMPID',all.y=T)
  outWide <- outMet.3

  # If we re-melt df now, we have missing values where the metric should be a
  # zero, so we can set NAs to 0 now
  outWide[is.na(outWide)] <- 0
  # outLong.1 <- data.table::melt(outMet.1,id.vars=c(sampID,'SAMPID')) %>%
  #   plyr::mutate(value=ifelse(is.na(value),0,value))

  # # Finally, we can recast the metrics df into wide format for output
  # lside <- paste(paste(sampID,collapse='+'),'SAMPID',sep='+')
  # formula <- paste(lside,'~variable',sep='')
  # outWide <- data.table::dcast(outLong.1,eval(formula),value.var='value')

  outAll <- merge(outWide,totals,by='SAMPID',all.x=T)
  outAll$SAMPID <- NULL
  # outAll <- merge(outWide,totals,by='SAMPID',all.x=T) %>%
  #   dplyr::select(-SAMPID)

  return(outAll)

}
