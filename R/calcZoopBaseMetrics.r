#' @export
#' @title Calculate base zooplankton metrics for given input dataset.
#'
#' @description This function calculates all base metrics associated
#' with a given input data frame using the inputs provided. These
#' inputs include the names of count, biomass, and density variables,
#' as well as a variable indicating distinctness of each taxon.
#'
#' @param indata A data frame containing, at minimum, the variables
#' specified in the arguments for sampID, dist, ct, and taxa_id,
#' @param sampID A character vector containing the names of all
#' variables in indf that specify a unique sample. If not specified,
#' the default is \emph{UID}
#' @param is_distinct A string with the name of the distinctness variable,
#' which is assumed to have only values of 0 or 1. If not specified,
#' the default is \emph{IS_DISTINCT}.
#' @param ct A string with the name of the count variable. This
#' variable should be numeric or able to be converted to numeric.
#' @param biomass A string with the name of the biomass variable. This
#' variable should be numeric or able to be converted to numeric.
#' @param density A string with the name of the density variable. This
#' variable should be numeric or able to be converted to numeric.
#' @param inTaxa a data frame containing taxonomic information,
#' including the variables PHYLUM, CLASS, SUBCLASS, ORDER, SUBORDER,
#' FAMILY, following the same taxonomy as the default taxa list,
#' \emph{zoopTaxa}. There should also be autecology traits with names
#' that match those for the arguments \emph{ffg}, \emph{clad_size},
#' and \emph{net_size}.
#' @param taxa_id A string with the name of the taxon ID variable
#' in \emph{inCts} that matches that in \emph{inTaxa}. The default
#' value is \emph{TAXA_ID}.
#' @param ffg A string with the name of the functional feeding group
#' variable in \emph{inTaxa}. Values used
#' in calculations include FILT, HERB, OMNI, PARA, PRED, representing
#' filterer, herbivore, omnivore, parasite, and predator,
#' respectively.
#' @param clad_size A string with the name of the variable in
#' \emph{inTaxa} indicating the size class of the cladoceran
#' taxon. Valid values are LARGE and SMALL.
#' @param net_size A string with the name of the variable in
#' \emph{inTaxa} indicating the net size class of a taxon.
#' Valid values are COARSE and FINE.
#' @param nativeMetrics A logical argument. TRUE indicates that
#' the subset of metrics based on native status should be
#' calculated. FALSE indicates that the full set of metrics
#' should be calculated. The default value is FALSE. If value
#' is TRUE, the input data should already be subset to
#' native taxa.
#' @return A data frame containing the variables in sampID and
#' the zooplankton metrics as additional variables. The output
#' metrics include number of individuals, percent individuals,
#' biomass, percent biomass, density, and percent density for
#' each of several groups based on traits \emph{inTaxa}. If
#' \emph{nativeMetrics} = TRUE, NAT is appended to metric names.

calcZoopBaseMetrics <- function(indata, sampID, is_distinct,
                           ct, biomass, density = NULL,
                           inTaxa, taxa_id='TAXA_ID',
                           ffg, clad_size, net_size,
                           nativeMetrics = FALSE){

  indata <- as.data.frame(indata)

  necVars <- c(sampID, is_distinct, ct, biomass, taxa_id)
  if(any(necVars %nin% names(indata))){
    msgTraits <- which(necVars %nin% names(indata))
    print(paste("Missing variables in input data frame:",
                paste(necVars[msgTraits], collapse=',')))
    return(NULL)
  }else{
    indata[, c(ct, biomass, is_distinct)] <- lapply(indata[, c(ct, biomass, is_distinct)], as.numeric)
  }

  if(!is.null(density)){
    if(density %nin% names(indata)){
      print("Missing density variable in input data frame.")
    }else{
      indata[, density] <- as.numeric(indata[, density])
    }
  }

  necTaxVars <- c(taxa_id, ffg, clad_size, net_size,
                  'PHYLUM', 'CLASS', 'SUBCLASS', 'ORDER',
                  'SUBORDER', 'FAMILY')
  if(any(necTaxVars %nin% names(inTaxa))){
    msgTraits <- which(necTaxVars %nin% names(inTaxa))
    print(paste("Missing variables in input taxalist:",
                paste(necTaxVars[msgTraits], collapse=',')))
    return(NULL)
  }

  # Assign characteristics to taxa
  calcData <- merge(indata, inTaxa, by = 'TAXA_ID')

  calcData$CALAN <- with(calcData, ifelse(ORDER=='CALANOIDA', 1, 0))
  calcData$COPE <- with(calcData, ifelse(SUBCLASS=='COPEPODA', 1, 0))
  calcData$COARSE <- with(calcData, ifelse(eval(as.name(net_size))=='COARSE', 1, 0))
  calcData$FINE <- with(calcData, ifelse(eval(as.name(net_size))=='FINE', 1, 0))
  calcData$SMCLAD <- with(calcData, ifelse(SUBORDER=='CLADOCERA' &
                                             eval(as.name(clad_size))=='SMALL', 1, 0))
  calcData$DAPHNIID <- with(calcData, ifelse(FAMILY=='DAPHNIIDAE', 1, 0))
  calcData$LGCLAD <- with(calcData, ifelse(eval(as.name(clad_size))=='LARGE', 1, 0))
  calcData$BOSM <- with(calcData, ifelse(FAMILY=='BOSMINIDAE', 1, 0))
  calcData$CLAD <- with(calcData, ifelse(SUBORDER=='CLADOCERA', 1, 0))

  # If input is not just native taxa, there are many more metrics that
  # get calculated.
  if(nativeMetrics==FALSE){
    calcData$COPE_HERB <- with(calcData, ifelse(COPE==1 & eval(as.name(ffg))=='HERB', 1, 0))
    calcData$ROT <- with(calcData, ifelse(PHYLUM=='ROTIFERA', 1, 0))
    calcData$HERB <- with(calcData, ifelse(eval(as.name(ffg))=='HERB', 1, 0))
    calcData$OMNI <- with(calcData, ifelse(eval(as.name(ffg))=='OMNI', 1, 0))
    calcData$PLOIMA <- with(calcData, ifelse(ORDER=='PLOIMA', 1, 0))
    calcData$SIDID <- with(calcData, ifelse(FAMILY=='SIDIDAE', 1, 0))
    calcData$CRUST <- with(calcData, ifelse(PHYLUM=='ARTHROPODA' &
                            CLASS %in% c('MAXILLOPODA','BRANCHIOPODA'), 1, 0))
    calcData$CYCLOP <- with(calcData, ifelse(ORDER=='CYCLOPOIDA', 1, 0))
    calcData$FLOS <- with(calcData, ifelse(ORDER=='FLOSCULARIACEAE', 1, 0))
    calcData$COLLO <- with(calcData, ifelse(ORDER=='COLLOTHECACEAE', 1, 0))
    calcData$ASPLAN <- with(calcData, ifelse(FAMILY=='ASPLANCHNIDAE', 1, 0))
    calcData$PRED <- with(calcData, ifelse(eval(as.name(ffg))=='PRED', 1, 0))
    calcData$CRUST_HERB <- with(calcData, ifelse(HERB==1 & CRUST==1, 1, 0))
    calcData$CRUST_PRED <- with(calcData, ifelse(PRED==1 & CRUST==1, 1, 0))
    calcData$CRUST_OMNI <- with(calcData, ifelse(OMNI==1 & CRUST==1, 1, 0))
    calcData$ROT_HERB <- with(calcData, ifelse(HERB==1 & ROT==1, 1, 0))
    calcData$ROT_PRED <- with(calcData, ifelse(PRED==1 & ROT==1, 1, 0))
    calcData$ROT_OMNI <- with(calcData, ifelse(OMNI==1 & ROT==1, 1, 0))
    calcData$COPE_PRED <- with(calcData, ifelse(PRED==1 & COPE==1, 1, 0))
    calcData$COPE_OMNI <- with(calcData, ifelse(OMNI==1 & COPE==1, 1, 0))
    calcData$CLAD_PRED <- with(calcData, ifelse(PRED==1 & CLAD==1, 1, 0))
    calcData$CLAD_OMNI <- with(calcData, ifelse(OMNI==1 & CLAD==1, 1, 0))
    calcData$CLAD_HERB <- with(calcData, ifelse(HERB==1 & CLAD==1, 1, 0))
  }

  # Now we create the list of parameters to run, and these depend on
  # whether nativeMetrics is TRUE or FALSE
  if(nativeMetrics == FALSE){
    params <- c('CALAN','COPE','COPE_HERB','ROT','COARSE','FINE','SMCLAD',
                'DAPHNIID','HERB','LGCLAD','OMNI','PLOIMA','SIDID','CLAD','BOSM'
                ,'CYCLOP','FLOS','COLLO','ASPLAN','PRED','ROT_HERB','ROT_PRED',
                'ROT_OMNI','COPE_PRED','COPE_OMNI','CLAD_HERB','CLAD_OMNI','CLAD_PRED')
  }else{
    params <- c('CALAN','COPE','COARSE','FINE','SMCLAD','LGCLAD','CLAD','DAPHNIID','BOSM')
  }

  samps <- unique(calcData[, sampID])

  column_names <- c(sampID, 'PARAMETER', 'RESULT')

  metsOut <- data.frame(matrix(nrow = 0, ncol = length(column_names)))
  colnames(metsOut) <- column_names

  # Need to calculate totals first
  if(!is.null(density)){
    totals <- calcZoopTotals(calcData, sampID, is_distinct,
                             c(ct, biomass, density),
                             c('TOTL_NIND', 'TOTL_BIO', 'TOTL_DEN'),
                             'TOTL_NTAX')
  }else{
    totals <- calcZoopTotals(calcData, sampID, is_distinct,
                             c(ct, biomass, density),
                             outputSums = c('TOTL_NIND', 'TOTL_BIO'),
                             outputTaxa = 'TOTL_NTAX')
  }


  for(i in 1:length(params)){
      print(i)
      metsIn <- subset(calcData, eval(as.name(params[i]))==1)

      metsIn <- merge(metsIn, totals, by = sampID)

      # For parameters with no observations, print a message for user
      # Then create an set of metrics for that parameter and assign value of 0
      # Then it adds onto the existing metrics calculated
      if(nrow(metsIn)==0){
        print(paste("No observations for ", params[i], sep=''))
        numSamps <- nrow(metsOut)
        met.1.long <- merge(samps, data.frame(PARAMETER=c(paste(params[i],'DEN', sep='_'),
                                               paste(params[i],'PDEN', sep='_'),
                                               paste(params[i],'BIO',sep='_'),
                                               paste(params[i],'PBIO',sep='_'),
                                               paste(params[i],'PIND',sep='_'),
                                               paste(params[i],'PTAX',sep='_'),
                                               paste(params[i],'NTAX',sep='_')),
                                              RESULT=0))
        metsOut <- rbind(metsOut, met.1.long)

      }else{ # Include NIND only for full metrics, not native only
        met.1 <- metsIn
        # If there are samples with parameter present, calculate all types of metrics for
        # that parameter, then add to existing metrics
        met.1$CALCPIND <- with(met.1, eval(as.name(ct))/TOTL_NIND)
        met.1$CALCPTAX <- with(met.1, eval(as.name(is_distinct))/TOTL_NTAX)
        met.1$CALCPDEN <- with(met.1, eval(as.name(density))/TOTL_DEN)
        met.1$CALCPBIO <- with(met.1, eval(as.name(biomass))/TOTL_BIO)

        met.1a.nind <- aggregate(x = list(NIND = met.1[, ct],
                                          NTAX = met.1[, is_distinct]),
                            by = met.1[c(sampID)],
                            FUN = function(x){sum(x, na.rm=T)})

        if(!is.null(density)){
          met.1a.den <- aggregate(x = list(DEN = met.1[, density]),
                                  by = met.1[c(sampID)],
                                  FUN = function(x){round(sum(x, na.rm=T), 4)})
        }

        met.1a.bio <- aggregate(x = list(BIO = met.1[, biomass]),
                                by = met.1[c(sampID)],
                                FUN = function(x){round(sum(x, na.rm=T), 6)})

        met.1b <- aggregate(x = list(PIND = met.1$CALCPIND,
                                     PTAX = met.1$CALCPTAX,
                                     PDEN = met.1$CALCPDEN,
                                     PBIO = met.1$CALCPBIO),
                            by = met.1[c(sampID)],
                            FUN = function(x){round(sum(x, na.rm=TRUE)*100, 2)})
        if(!is.null(density)){
          met.2 <- merge(met.1a.nind, met.1a.bio, by = sampID) |>
            merge(met.1a.den, by = sampID) |>
            merge(met.1b, by = sampID)
        }else{
          met.2 <- merge(met.1a.nind, met.1a.bio, by = sampID) |>
            merge(met.1b, by = sampID)
        }

        met.2a <- merge(samps, met.2, by = c(sampID), all.x=TRUE)
        met.2a[is.na(met.2a)] <- 0

        met.2.long <- reshape(met.2a, idvar = sampID, direction = 'long',
                            varying = names(met.2)[!(names(met.2a) %in% sampID)],
                            timevar = 'PARAMETER', v.names = 'RESULT',
                            times = names(met.2)[!(names(met.2a) %in% sampID)])

        if(nativeMetrics==FALSE){
          met.2.long <- mutate(met.2.long, PARAMETER=paste(params[i], PARAMETER, sep='_'),
                             RESULT=ifelse(is.na(RESULT), 0, RESULT))
        }else{
          met.2.long <- mutate(met.2.long, PARAMETER=paste(params[i], 'NAT', PARAMETER, sep='_'),
                               RESULT=ifelse(is.na(RESULT), 0, RESULT))
        }

        metsOut <- rbind(metsOut, met.2.long)

      }
  }
  return(metsOut)
}

