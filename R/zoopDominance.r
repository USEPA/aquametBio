#' @export
#' @keywords internal
#' @title Dominance metric calculation function (zooplankton)
#' @description This function calculates % dominant organisms metric
#' for zooplankton. This is done slightly differently from other
#' assemblages, including rounding to 1 digit instead of 2. It also
#' accepts different input variables, which accommodates calculating
#' dominance based on count, biomass, and density.
#' @param df Input data frame, containing SAMPID as variable identifying unique samples
#' @param sampID A character vector containing the names of all
#' variables in indf that specify a unique sample. If not specified,
#' the default is \emph{UID}.
#' @param topN Number specifying the top number of species to include in calculation.
#' The default value is 1.
#' @param varIn A string with the name of the variable to be used in calculating
#' dominance.
#' @return A data frame with \emph{sampID} variables and the metric containing
#' the % individuals in the
#' dominant (topN) taxa
zoopDominance<-function(df, sampID = 'UID', topN = 1, varIn){

  df.long <- reshape(df, idvar = sampID, direction = 'long',
                     varying = varIn, timevar = 'variable',
                     v.names = 'value', times = varIn)

  rr <- aggregate(list(TOTSUM = df.long$value),
                  list(df.long[, sampID]),
                  function(x){sum(x, na.rm=TRUE)})

  ss <- merge(df.long, rr, by = sampID)

  tt <- aggregate(list(domN = ss$value),
                  list(ss[, sampID]),
                  function(x){
                    sum(x[order(x, decreasing=TRUE)[1:topN]],
                        na.rm=T)
                  }
  )
  uu <- merge(tt, unique(ss[, c(sampID,'TOTSUM')]), by = sampID)
  uu <- mutate(uu, dompind=round(domN/TOTSUM*100, 1))
  uu <- subset(uu, select=c(sampID, 'dompind'))

  return(uu)
}
