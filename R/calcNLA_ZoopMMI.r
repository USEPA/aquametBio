#' @export
#' @title Calculate the NLA zooplankton MMI
#'
#' @description This is a function that calculates
#' the benthic MMI as used for the National Lakes
#' Assessment, based on inputs of the appropriate metrics.
#'
#' @param inMets A data frame containing, at minimum, the variables
#' specified in the arguments for sampID, ecoreg, totlnind, and the metrics
#' necessary for calculation of MMI by region. If values for more than
#' the necessary metrics are included in the input data frame, the
#' unnecessary  metrics will be ignored for each given site.
#'
#' The necessary metrics, by aggregated bioregion, are:
#'
#'  CPL: FINE_BIO, SIDID_PIND, DOM1_300_COPE_PBIO, FAM300_NAT_NTAX,
#'  COLLO_PBIO, OMNI_PTAX
#'
#'  EHIGH: ZOCN_DEN, SMCLAD_PBIO, COPE_NAT_DEN, COARSE_NAT_PTAX, ROT_PBIO,
#'  OMNI300_PTAX
#'
#'  PLAINS: FINE300_NAT_PBIO, SMCLAD_NAT_PIND, COPE_RATIO_300_BIO,
#'  FAM300_NAT_NTAX, ROT_NTAX, COPE_HERB_PDEN
#'
#'  UMW: TOTL_NAT_PIND, BOSM300_NAT_PTAX, CALAN300_NAT_BIO, FINE_PTAX,
#'  DOM1_ROT_PBIO, COPE_HERB300_PBIO
#'
#'  WMTNS: COARSE300_NAT_PBIO, LGCLAD300_NAT_PTAX, COPE300_BIO,
#'  ZOFN300_NTAX, PLOIMA_PTAX, COPE_OMNI_PTAX
#'
#'  Descriptions of these metrics can be found in the file
#'  \emph{NRSA_Invertebrate_Metric_Descriptions.pdf}, included in
#'  the documentation for this package.
#'
#' @param sampID A character vector containing the names of all
#' variables in inMets that specify a unique sample. If not specified,
#' the default is \emph{UID}
#' @param ecoreg A string with the name of the aggregated bioregion variable.
#' Valid values that correspond to regions used in NLA are
#' CPL, EHIGH, PLAINS, UMW, and WMTNS.
#' @param totlnind A string with the name of the variable with the
#' total individuals in each sample.
#' @return A data frame containing the variables in sampID, as well as
#' the scored metrics, the zooplankton MMI, and the condition class for each
#' sites. The variable names are ABUN_PT, CLAD_PT, COPE_PT, RICH_PT,
#' ROTI_PT, TROP_PT, MMI_ZOOP, ZOOP_MMI_COND.
#' @author Karen Blocksom \email{Blocksom.Karen@epa.gov}
#' @keywords survey
calcNLA_ZoopMMI <- function(inMets, sampID='UID', ecoreg='ECOREG', totlnind = 'TOTL_NIND'){
  # Convert data into data frames just in case
  inMets <- as.data.frame(inMets)

  necTraits <- c(sampID, ecoreg, totlnind)
  if(any(necTraits %nin% names(inMets))){
    msgTraits <- which(necTraits %nin% names(inMets))
    print(paste("Some of the traits are missing from the input dataset. The following are required for metric calculations to run:"
                , necTraits[msgTraits]))
    return(NULL)
  }

  # Rename variables
  names(inMets)[names(inMets)==ecoreg] <- 'ECO_BIO'
  names(inMets)[names(inMets)==totlnind] <- 'TOTLNIND'

  # Combine all values in sampID into one sampID in df
  for(i in 1:length(sampID)){
    if(i==1) inMets$SAMPID <- inMets[,sampID[i]]
    else inMets$SAMPID <- paste(inMets$SAMPID,inMets[,sampID[i]],sep='.')
  }
  samples <- unique(inMets[,c('SAMPID', sampID)])

  # Check to make sure ecoregion variable is included in the input data frame
  ecoCk <- unique(inMets$ECO_BIO)
  ecos <- c('CPL','EHIGH','PLAINS','UMW','WMTNS')
  if(any(ecoCk %nin% ecos)){
    msgEco <- which(ecoCk %nin% ecos)
    print(paste("These ecoregions are not valid: "
                ,paste(ecoCk[msgEco],collapse=',')))
    return(NULL)
  }

  cfVal <- data.frame(ECO_BIO=c(rep('CPL',6), rep('EHIGH',6), rep('PLAINS',6), rep('UMW',6), rep('WMTNS',6)),
                      PARAMETER=c('FINE_BIO','SIDID_PIND','DOM1_300_COPE_PBIO','FAM300_NAT_NTAX',
                                  'COLLO_PBIO','OMNI_PTAX',
                                  'ZOCN_DEN','SMCLAD_PBIO','COPE_NAT_DEN',
                                  'COARSE_NAT_PTAX','ROT_PBIO','OMNI300_PTAX',
                                  'FINE300_NAT_PBIO','SMCLAD_NAT_PIND','COPE_RATIO_300_BIO',
                                  'FAM300_NAT_NTAX','ROT_NTAX','COPE_HERB_PDEN',
                                  'TOTL_NAT_PIND','BOSM300_NAT_PTAX','CALAN300_NAT_BIO',
                                  'FINE_PTAX','DOM1_ROT_PBIO','COPE_HERB300_PBIO',
                                  'COARSE300_NAT_PBIO','LGCLAD300_NAT_PTAX','COPE300_BIO',
                                  'ZOFN300_NTAX','PLOIMA_PTAX','COPE_OMNI_PTAX'),
                      MET_TYPE=rep(c('ABUN','CLAD','COPE','RICH','ROTI','TROP'),5),
                      DISTRESP=c('NEGATIVE','NEGATIVE','POSITIVE','POSITIVE','POSITIVE','NEGATIVE',
                                 'NEGATIVE','NEGATIVE','NEGATIVE','POSITIVE','NEGATIVE','NEGATIVE',
                                 'POSITIVE','POSITIVE','POSITIVE','POSITIVE','POSITIVE','NEGATIVE',
                                 'POSITIVE','POSITIVE','NEGATIVE','POSITIVE','NEGATIVE','NEGATIVE',
                                 'POSITIVE','POSITIVE','NEGATIVE','NEGATIVE','POSITIVE','NEGATIVE'),
                      FLOOR=c(2.913623, 0, 45.9, 5, 0, 10.53,
                              0.21645, 0, 8.8236, 22.22, 0.79, 12.5,
                              0.66, 0, 0, 5, 3, 0,
                              96.75, 0, 0, 37.5, 25.3, 0.19,
                              10.94, 0, 0.07393, 3, 20, 0),
                      CEILING=c(173.28, 24.88, 100, 15, 5.64, 47.06,
                                259.305, 57.31, 398.397, 57.14, 86.39, 44.44,
                                85.12, 49.03, 62.8097, 15, 17, 29.93,
                                100, 12.5, 65.0375, 77.78, 93.6, 59.42,
                                99.26, 29.285, 149.036, 15, 70.835, 22.22))

  matchMets <- reshape(inMets, idvar = c('SAMPID','ECO_BIO','TOTLNIND'),
                       direction = 'long',
                       varying = names(inMets)[names(inMets) %in% unique(cfVal$PARAMETER)],
                       timevar = 'PARAMETER', v.names = 'RESULT',
                       times = names(inMets)[names(inMets) %in% unique(cfVal$PARAMETER)])

  matchMets <- merge(matchMets, cfVal, by = c('PARAMETER', 'ECO_BIO'))

  # Run a check to make sure there are exactly 6 rows per sites in the matched dataset
  numMets <- as.data.frame(table(SAMPID = matchMets$SAMPID))
  numMets <- subset(numMets, Freq<6)
  numMets <- merge(numMets, inMets, by = 'SAMPID')

  if(nrow(numMets)>0){
    return(print(paste0("Missing metrics values for these samples: ",
                       numMets$SAMPID,
                       ". Check input data frame against required metric list.")))
  }


  ## The function below interpolates the score between the floor and ceiling scoring thresholds for each metric
  scoreMet1 <- function(resptype, x, floor, ceiling){
    if(resptype=='POSITIVE'){
      zz<-round(approx(x=c(floor, ceiling), y=c(0, 10), xout=x, method='linear', yleft=0, yright=10)$y, 2)
    } else {
      zz<-round(approx(x=c(floor, ceiling), y=c(10, 0),xout=x, method='linear', yleft=10, yright=0)$y, 2)
    }

  }

  ## Send metric values to the scoring function above (scoreMet1)
  scored.mets <- matchMets[,c('SAMPID', 'TOTLNIND', 'ECO_BIO','MET_TYPE', 'PARAMETER')]
  scored.mets$RESULT <- with(scored.mets, ifelse(as.numeric(TOTLNIND)==0, 0,
                                      with(matchMets, mapply(scoreMet1, DISTRESP, RESULT, FLOOR, CEILING))))
  scored.mets$PARAMETER <- paste(scored.mets$MET_TYPE, 'PT', sep = '_')

  ## Sum metrics scores for each sample and rescale total to 100-point scale
  mmi.scores <- aggregate(x = list(SUMMETS = scored.mets$RESULT), by = scored.mets[c('SAMPID','TOTLNIND','ECO_BIO')], FUN = sum)

  mmi.scores$PARAMETER <- "MMI_ZOOP"
  mmi.scores$RESULT <- with(mmi.scores, round((100/60)*SUMMETS, 1))
  mmi.scores$SUMMETS <- NULL

  ## Set condition class for each sample, which is based on ECO_BIO region
  # First create a table of thresholds by ECO_BIO
  condTholds <- data.frame(ECO_BIO = c('EHIGH', 'CPL', 'WMTNS', 'PLAINS', 'UMW'),
                           gf = c(73.595, 59.42, 60.78, 36.72, 63.68),
                           fp = c(60.03, 53.77, 51.32, 28.17, 52.03))

  ## Merge MMI scores with thresholds by ECO9 region
  cond.mmi <- merge(mmi.scores,condTholds, by='ECO_BIO')
  cond.mmi$ECO9 <- as.character(cond.mmi$ECO_BIO)
  cond.mmi$PARAMETER <- 'ZOOP_MMI_COND'
  cond.mmi$MMI_ZOOP <- cond.mmi$RESULT
  cond.mmi$RESULT <- with(cond.mmi, ifelse(is.na(MMI_ZOOP), 'Not Assessed',
                                           ifelse(MMI_ZOOP>=gf, 'Good',
                                                  ifelse(MMI_ZOOP<fp, 'Poor', 'Fair'))))

  ww <- rbind(subset(scored.mets, select=c('SAMPID', 'ECO_BIO', 'PARAMETER', 'RESULT'))
              ,subset(mmi.scores, select=c('SAMPID', 'ECO_BIO', 'PARAMETER', 'RESULT'))
              ,subset(cond.mmi, select=c('SAMPID', 'ECO_BIO', 'PARAMETER', 'RESULT')))

  # Finally, we can recast the metrics df into wide format for output
  mmiOut <- reshape(ww, idvar = c('SAMPID', 'ECO_BIO'), direction = 'wide',
                    v.names = 'RESULT', timevar = 'PARAMETER')
  names(mmiOut) <- gsub('RESULT\\.', '', names(mmiOut))

  mmiOut.final <- merge(samples, mmiOut, by='SAMPID')
  mmiOut.final <- subset(mmiOut.final, select=c(sampID, 'SAMPID', 'ECO_BIO', 'MMI_ZOOP',
                                                'ZOOP_MMI_COND',
                                                names(mmiOut)[names(mmiOut) %nin% c(sampID, 'SAMPID',
                                                                                    'ECO_BIO',
                                                                                    'MMI_ZOOP',
                                                                                    'ZOOP_MMI_COND')]))
  names(mmiOut.final)[names(mmiOut.final)=='ECO_BIO'] <- ecoreg

  mmiOut.final$ZOOP_MMI_COND <- with(mmiOut.final, ifelse(is.na(MMI_ZOOP), 'Not Assessed', ZOOP_MMI_COND))

  return(mmiOut.final)


}
