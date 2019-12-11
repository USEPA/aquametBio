#' @export
#' @title Calculate anomaly-based metrics
#' @description This function calculates all of the metrics based
#' only on fish anomalies.
#' @param indata Input data frame containing variables as identified
#' in the arguments for \emph{sampID}, \emph{dist}, \emph{ct},
#' \emph{taxa_id}.
#' @param sampID sampID A character vector containing the names of all
#' variables in indf that specify a unique sample. If not specified,
#' the default is \emph{UID}.
#' #' @param ct A string with the name of the count variable. If not
#' specified, the default is \emph{TOTAL}.
#' @param ct A string with the name of the count variable. If not
#' specified, the default is \emph{TOTAL}.
#' @param anomct A string with the name of the count variable. If not
#' specified, the default is \emph{ANOM_CT}.
#' @return A data frame containing the variables in sampID and
#' the fish anomaly metric ANOMPIND as an additional variable. Metric
#' descriptions are included in \emph{NRSA_Fish_Metric_Descriptions.pdf},
#' included in this package.
#'
#' @author Karen Blocksom \email{Blocksom.Karen@epa.gov}
#'
calcFishAnomMets <- function(indata, sampID='UID',
                             ct='TOTAL', anomct='ANOM_CT'){

  ctVars <- c(sampID,ct,anomct)
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
  names(indata)[names(indata)==anomct] <- 'ANOM_CT'

  indata <- plyr::mutate(indata,ANOM_CT=as.numeric(ANOM_CT),TOTAL=as.numeric(TOTAL))

  # Make sure ANOM_CT value not present when TOTAL is missing or 0
  checkAnom <- subset(indata,is.na(TOTAL)|TOTAL==0 & !is.na(ANOM_CT) & ANOM_CT!=0)
  if(nrow(checkAnom)>0){
    print(head(checkAnom))
    return(print("Values for anomalies non-zero when count for species missing or zero. Cannot calculate metric correctly."))
  }

  outMet <- aggregate(x = list(ANOMNIND = indata$ANOM_CT, TOTLNIND = indata$TOTAL), by = indata[c('SAMPID')], FUN = function(x){sum(x, na.rm=T)})
  outMet$ANOMPIND <- with(outMet, round(ANOMNIND/TOTLNIND*100, 2))
  outMet <- subset(outMet, select = c('SAMPID','ANOMPIND'))
  # outMet <- plyr::ddply(indata, "SAMPID", summarise, ANOMPIND=round(sum(ANOM_CT,na.rm=T)/sum(TOTAL,na.rm=T)*100,2))

  outMet.1 <- merge(samples, outMet, by = 'SAMPID', all.x=TRUE)
  outMet.1$ANOMPIND <- with(outMet.1, ifelse(is.na(ANOMPIND),0,ANOMPIND))
  # outMet.1 <- merge(samples,outMet,by='SAMPID') %>%
  #   plyr::mutate(ANOMPIND=ifelse(is.na(ANOMPIND),0,ANOMPIND))

  outAll <- subset(outMet.1, select = -SAMPID)
  # outAll <- dplyr::select(outMet.1,-SAMPID)

  return(outAll)

}
