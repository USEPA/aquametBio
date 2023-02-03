#' @export
#' @title Calculate zooplankton diversity metrics for given input dataset.
#'
#' @description This function calculates all diversity metrics associated
#' with a given input data frame using the inputs provided.
#'
#' @param indata A data frame containing, at minimum, the variables
#' specified in the arguments for sampID, dist, ct, and taxa_id
#' and TRIBE, as well as autecology traits with names that match those
#' in the arguments ffg, habit, and ptv. In addition, there
#' should be a variable with the name in argument taxa_id that matches
#' with all of those in the indf data frame
#' @param sampID A character vector containing the names of all
#' variables in indata that specify a unique sample. If not specified,
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
#' @param suffix A string to indicate the suffix that should be added
#' to metric names to indicate the subgroup represented in the metric.
#' This suffix must match the name of a variable in \emph{indata}
#' @return A data frame containing the variables in sampID and
#' the zooplankton metrics as additional variables. If
#' \emph{nativeMetrics} = TRUE, NAT is appended to metric names.
calcZoopDivMetrics <- function(indata, sampID, is_distinct,
                            ct, biomass = NULL, density = NULL,
                            suffix = ''){
  indata <- as.data.frame(indata)

  necVars <- c(sampID, is_distinct, ct)
  if(any(necVars %nin% names(indata))){
    msgTraits <- which(necVars %nin% names(indata))
    print(paste("Missing variables in input data frame:",
                paste(necVars[msgTraits], collapse=',')))
    return(NULL)
  }else{
    indata[, c(ct, is_distinct)] <- lapply(indata[, c(ct, is_distinct)], as.numeric)
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
  hprime <- aggregate(x=list(HPRIME_NIND = calcData[, 'prop.cnt']),
                             by = calcData[sampID],
                             FUN = function(x){round(-1*sum(x*log(x), na.rm=T), 4)})

  simpson <- aggregate(x = list(SIMPSON_NIND = calcData[, 'prop.cnt']),
                       by = calcData[sampID],
                       FUN = function(x){round(sum(x*x, na.rm=T), 4)})

  even <- merge(hprime, totals[, c(sampID, 'TOTL_NTAX')], by=sampID)
  even$EVEN_NIND <- with(even, round(HPRIME_NIND/log(TOTL_NTAX), 4))

  even.1 <- subset(even, select = c(sampID, 'EVEN_NIND'))

  pie.in <- calcData
  pie.in$CALCPROP <- with(pie.in, prop.cnt*((TOTL_NIND - eval(as.name(ct)))/(TOTL_NIND-1)))

  pie <- aggregate(x=list(PIE_NIND = pie.in[, 'CALCPROP']),
              by = pie.in[sampID],
              FUN = function(x){round(sum(x, na.rm=T), 4)})

  metsOut <- merge(hprime, simpson, by = sampID) |>
    merge(pie, by = sampID) %>%
    merge(even.1, by = sampID)

  if(!is.null(density)){
    hprime.den <- aggregate(x=list(HPRIME_DEN = calcData[, 'prop.den']),
                        by = calcData[sampID],
                        FUN = function(x){round(-1*sum(x*log(x), na.rm=T), 4)})

    simpson.den <- aggregate(x = list(SIMPSON_DEN = calcData[, 'prop.den']),
                         by = calcData[sampID],
                         FUN = function(x){round(sum(x*x, na.rm=T), 4)})

    metsOut <- merge(metsOut, hprime.den, by = sampID) |>
      merge(simpson.den, by = sampID)
  }

  if(!is.null(biomass)){
    hprime.bio <- aggregate(x=list(HPRIME_BIO = calcData[, 'prop.bio']),
                        by = calcData[sampID],
                        FUN = function(x){round(-1*sum(x*log(x), na.rm=T), 4)})

    simpson.bio <- aggregate(x = list(SIMPSON_BIO = calcData[, 'prop.bio']),
                         by = calcData[sampID],
                         FUN = function(x){round(sum(x*x, na.rm=T), 4)})

    metsOut <- merge(metsOut, hprime.bio, by = sampID) |>
      merge(simpson.bio, by = sampID)

  }

  metsOut.long <- reshape(metsOut, idvar = sampID, direction = 'long',
                        varying = names(metsOut)[!(names(metsOut) %in% sampID)],
                        timevar = 'PARAMETER', v.names = 'RESULT',
                        times = names(metsOut)[!(names(metsOut) %in% sampID)])

  metsOut.long$RESULT <- ifelse(is.na(metsOut.long$RESULT), 0, metsOut.long$RESULT)

  if(suffix != ''){
    metsOut.long$PARAMETER <- gsub('NIND', suffix, metsOut.long$PARAMETER)
  }

  return(metsOut.long)

}
