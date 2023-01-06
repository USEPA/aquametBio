#' @export
#' @title Calculate zooplankton diversity metrics for given input dataset.
#'
#' @description This function calculates all diversity metrics associated
#' with a given input data frame using the inputs provided.
#'
#' @param inCts A data frame containing, at minimum, the variables
#' specified in the arguments for sampID, dist, ct, and taxa_id
#' and TRIBE, as well as autecology traits with names that match those
#' in the arguments ffg, habit, and ptv. In addition, there
#' should be a variable with the name in argument taxa_id that matches
#' with all of those in the indf data frame
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
#' the zooplankton metrics as additional variables. If
#' \emph{nativeMetrics} = TRUE, NAT is appended to metric names.
calcZoopDivMetrics <- function(indata, sampID, is_distinct,
                            ct, biomass = NULL, density = NULL){
  indata <- as.data.frame(indata)

  necVars <- c(sampID, is_distinct, ct, taxa_id)
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

  if(!is.null(biomass)){
    if(biomass %nin% names(indata)){
      print("Missing biomass variable in input data frame.")
    }else{
      indata[, biomass] <- as.numeric(indata[, biomass])
    }
  }

  indata.1 <- subset(indata, eval(as.name(is_distinct))==1 & eval(as.name(ct))>0 &
                       !is.na(eval(as.name(ct))))
  ## Need to recalculate total taxa because TOTL_NTAX includes large rare and invasive taxa without counts associated
  totals <- calcZoopTotals(indata.1, sampID, is_distinct,
                              ct, 'TOTL_NIND', 'TOTL_NTAX')

  calcData <- merge(indata.1, totals, by = sampID)
  calcData$prop.cnt <- calcData[, ct]/calcData$TOTL_NIND

  if(!is.null(density)){
    totals.dens <- calcZoopTotals(indata.1, sampID, is_distinct,
                             c(density),
                             c('TOTL_DEN'))

    calcData <- merge(calcData, totals.dens, by = sampID)
    calcData$prop.den <- calcData[, density]/calcData$TOTL_DEN

  }

  if(!is.null(biomass)){
    totals.bio <- calcZoopTotals(indata.1, sampID, is_distinct,
                             c(biomass),
                             outputSums = c('TOTL_BIO'))

    calcData <- merge(calcData, totals.bio, by = sampID)
    calcData$prop.bio <- calcData[, biomass]/calcData$TOTL_BIO

  }

  # Calculate diversity for counts, biomass, and density separately to
  # accommodate all variations of calculation of indices
  hprime <- aggregate(x=list(HPRIME_NIND = calcData[, prop.cnt]),
                             by = calcData[, c(sampID)],
                             FUN = function(x){round(-1*sum(x*log(x), na.rm=T), 4)})

  simpson <- aggregate(x = list(SIMPSON_NIND = calcData[, prop.cnt]),
                       by = calcData[, c(sampID)],
                       FUN = function(x){round(sum(x*x, na.rm=T), 4)})

  even <- merge(hprime, totals[, c(sampID, 'TOTL_NTAX')], by=sampID)
  even$EVEN_NIND <- with(even, round(HPRIME_NIND/log(TOTL_NTAX), 4))

  pie.in <- calcData
  pie.in$CALCPROP <- with(pie.in, prop.cnt*((TOTL_NIND - COUNT)/(TOTL_NIND-1)))

  pie <- aggregate(x=list(PIE_NIND = pie.in[, 'CALCPROP']),
              by = pie.in[, c(sampID)],
              FUN = function(x){round(sum(x, na.rm=T), 4)})

  metsOut <- merge(hprime, simpson, by = sampID) |>
    merge(pie, by = sampID)

  if(!is.null(density)){
    hprime.den <- aggregate(x=list(HPRIME_DEN = calcData[, prop.den]),
                        by = calcData[, c(sampID)],
                        FUN = function(x){round(-1*sum(x*log(x), na.rm=T), 4)})

    simpson.den <- aggregate(x = list(SIMPSON_DEN = calcData[, prop.den]),
                         by = calcData[, c(sampID)],
                         FUN = function(x){round(sum(x*x, na.rm=T), 4)})

    metsOut <- merge(metsOut, hprime.den, by = sampID) |>
      merge(simpson.den, by = sampID)
  }

  if(!is.null(biomass)){
    hprime.bio <- aggregate(x=list(HPRIME_BIO = calcData[, prop.bio]),
                        by = calcData[, c(sampID)],
                        FUN = function(x){round(-1*sum(x*log(x), na.rm=T), 4)})

    simpson.bio <- aggregate(x = list(SIMPSON_BIO = calcData[, prop.bio]),
                         by = calcData[, c(sampID)],
                         FUN = function(x){round(sum(x*x, na.rm=T), 4)})

    metsOut <- merge(metsOut, hprime.bio, by = sampID) |>
      merge(simpson.bio, by = sampID)

  }


  # div <- ddply(indf.props, c('UID'), summarise,
  #              HPRIME_NIND=round(-1*sum(prop.cnt*log(prop.cnt),na.rm=T),4),
  #              HPRIME_BIO=round(-1*sum(prop.bio*log(prop.bio),na.rm=T),4),
  #              HPRIME_DEN=round(-1*sum(prop.den*log(prop.den),na.rm=T),4),
  #              SIMPSON_NIND=round(sum(prop.cnt*prop.cnt,na.rm=T),4),
  #              SIMPSON_BIO=round(sum(prop.bio*prop.bio,na.rm=T),4),
  #              SIMPSON_DEN=round(sum(prop.den*prop.den,na.rm=T),4),
  #              PIE_NIND=round(sum(prop.cnt*((SUMCNT-COUNT)/(SUMCNT-1)),na.rm=T),4),
  #              EVEN_NIND=round(HPRIME_NIND/log(unique(SUBTOTL_NTAX)),4))

}
