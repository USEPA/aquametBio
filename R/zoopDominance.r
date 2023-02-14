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
#' @param taxa_id A string with the name of the variable that distinctly
#' identifies taxa in each sample.
#' @return A data frame with \emph{sampID} variables and the metric containing
#' the % individuals in the
#' dominant (topN) taxa
zoopDominance <- function(df, sampID = 'UID', topN = 1, varIn, taxa_id){

  df.long <- reshape(df, idvar = c(sampID, taxa_id), direction = 'long',
                     varying = varIn, timevar = 'variable',
                     v.names = 'value', times = varIn) |>
    subset(select = c(sampID, taxa_id, 'variable', 'value'))

  rr <- aggregate(x = list(TOTSUM = df.long[, 'value']),
                  by = df.long[sampID],
                  FUN = function(x){sum(x, na.rm=TRUE)})

  ss <- merge(df.long, rr, by = sampID)

  tt <- aggregate(x = list(domN = ss[, 'value']),
                  by = ss[sampID],
                  function(x){
                    sum(x[order(x, decreasing=TRUE)[1:topN]],
                        na.rm=T)
                  }
  )
  uu <- merge(tt, unique(ss[, c(sampID,'TOTSUM')]), by = sampID)
  uu$dompind <- with(uu, round(domN/TOTSUM*100, 1))
  uu <- subset(uu, select=c(sampID, 'dompind'))

  return(uu)
}
