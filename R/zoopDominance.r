#' @export
#' @keywords internal
#' @title Dominance metric calculation function (zooplankton)
#' @description This function calculates % dominant organisms metric
#' for zooplankton. This is done slightly differently from other
#' assemblages, including rounding to 1 digit instead of 2. It also
#' accepts different input variables, which accommodates calculating
#' dominance based on count, biomass, and density.
#' @param df Input data frame, containing SAMPID as variable identifying unique samples
#' @param topN Number specifying the top number of species to include in calculation
#' @return A data frame with \emph{sampID} variables and the metric containing
#' the % individuals in the
#' dominant (topN) taxa
zoopDominance<-function(df, topN = 1, varIn){
  rr <- subset(df, IS_DISTINCT == 1)
  rr.long <- reshape(rr, )
  rr.long <- melt(rr, id.vars='UID' , measure.vars=varIn)

  ss <- ddply(rr.long, "UID", mutate, TOTSUM=sum(value, na.rm=T))

  tt <- aggregate(list(domN=ss$value)
                  ,list(UID=ss$UID)
                  ,function(x){
                    sum(x[order(x, decreasing=TRUE)[1:topN]],
                        na.rm=T)
                  }
  )
  uu <- merge(tt, unique(ss[, c('UID','TOTSUM')]), by="UID")
  uu <- mutate(uu, dompind=round(domN/TOTSUM*100, 1))
  uu <- subset(uu, select=c('UID', 'dompind'))

  return(uu)
}
