#' @export
#' @title Calculate zooplankton copepod ratio metrics
#' @description This function calculates ratios of calanoid
#' copepods to the sum of cyclopoid copepods and cladocerans
#' @param df Input data frame, containing SAMPID as variable identifying unique samples
#' @param sampID A character vector containing the names of all
#' variables in indf that specify a unique sample. If not specified,
#' the default is \emph{UID}.
#' @param is_distinct A string with the name of the distinctness variable,
#' which is assumed to have only values of 0 or 1. If not specified,
#' the default is \emph{IS_DISTINCT}.
#' @param valsIn A vector with the names of numeric variables to be used
#' in calculating dominance.
#' @param valsOut A vector of output suffixes to append to dominance
#' metric, matched with the order of valsIn variables.
#' @param taxa_id A string with the name of the variable that distinctly
#' identifies taxa in each sample.
#' @param subgrp A string with the name of a subgroup in \emph{indata}, with
#' a valid value of 1 to indicate inclusion in the subgroup.
#' @return A data frame with \emph{sampID} variables and the metric containing
#' the % individuals in the
#' dominant (topN) taxa
calcZoopCopeMetrics <- function(indata, sampID, is_distinct,
                               valsIn, valsOut, taxa_id,
                               subgrp = NULL){
}
