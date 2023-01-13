#' @export
#' @title Calculate percent native metrics
#' @description This function calculates percent native for
#' each of the variables in \emph{inputVars}.
#' @param indata Input data frame containing variables as identified
#' in the arguments for \emph{sampID} and \emph{inputVars}.
#' @param sampID sampID A character vector containing the names of all
#' variables in indf that specify a unique sample. If not specified,
#' the default is \emph{UID}.
#' @param is_distinct A string with the name of the distinctness variable,
#' which can have valid values of 0 (non-distinct) or 1 (distinct).
#' @param nonnative A string with the name of the numeric variable
#' indicating a non-native taxon. A value of 1 indicates that a
#' taxon is non-native. All other values will be ignored.
#' @param intaxa A string with the name of the data frame containing
#' the taxonomic data.
#'
#' @return A data frame containing percent native metrics
#' for the full data and the 300-count subsample
#'
#' @author Karen Blocksom \email{Blocksom.Karen@epa.gov}
#'
calcZoopNativeMetrics <- function(indata, sampID, is_distinct,
                                    inputVars, nonnative){

}
