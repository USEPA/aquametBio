#' @export
#' @title Calculate zooplankton copepod ratio metrics
#' @description This function calculates ratios of calanoid
#' copepods to the sum of cyclopoid copepods and cladocerans
#' @param indata Input data frame, containing \emph{sampID}
#' as variable identifying unique samples, as well as the
#' metrics listed in \emph{calaIn}, \emph{cyclIn}, and
#' \emph{cladIn}.
#' @param sampID A character vector containing the names of all
#' variables in \emph{indata} that specify a unique sample.
#' @param calaIn A vector of the names of input Calanoida metrics, in
#' the same order as for cladIn and cyclIn metrics.
#' @param cyclIn A vector of the names of input Cyclopoida metrics, in
#' the same order as for calaIn and cladIn metrics.
#' @param cladIn A vector of the names of input Cladocera metrics, in
#' the same order as for calaIn and cyclIn metrics.
#' @return A data frame with \emph{sampID} variables and the metric
#' containing the copepod ratio metrics
#'
#' \href{https://github.com/USEPA/aquametBio/blob/main/inst/NLA_Zooplankton_Metric_Descriptions.pdf}{NLA_Zooplankton_Metric_Descriptions.pdf}
#'
#' @author Karen Blocksom \email{Blocksom.Karen@epa.gov}
calcZoopCopeMetrics <- function(indata, sampID,
                                calaIn, cyclIn, cladIn) {
  indata <- as.data.frame(indata)

  if (length(calaIn) != length(cyclIn) |
    length(calaIn) != length(cladIn)) {
    print("calaIn, cyclIn, and cladIn must all contain the same
          number of variables.")
    return(NULL)
  }

  necVars <- c(sampID, calaIn, cyclIn, cladIn)
  if (any(necVars %nin% names(indata))) {
    msgTraits <- which(necVars %nin% names(indata))
    print(paste(
      "Missing variables in input data frame:",
      paste(necVars[msgTraits], collapse = ",")
    ))
    return(NULL)
  } else {
    indata[, c(calaIn, cyclIn, cladIn)] <- lapply(indata[, c(calaIn, cyclIn, cladIn)], as.numeric)
  }
  # Fill in zeros where missing values exist - this follows
  # what was done in NLA 2017
  indata[, c(calaIn, cyclIn, cladIn)] <- lapply(indata[, c(calaIn, cyclIn, cladIn)],
    FUN = function(x) {
      ifelse(is.na(x), 0, x)
    }
  )

  samps <- unique(indata[, sampID])

  column_names <- c(sampID, "PARAMETER", "RESULT")

  metsOut <- data.frame(matrix(nrow = 0, ncol = length(column_names)))
  colnames(metsOut) <- column_names

  for (i in 1:length(calaIn)) {
    temp <- indata
    temp$RESULT <- temp[, calaIn[i]] / (temp[, cyclIn[i]] + temp[, cladIn[i]])
    temp$RESULT <- ifelse(is.na(temp$RESULT) | is.infinite(temp$RESULT), 0, temp$RESULT)
    # to rename each one by the suffix on the input variables,
    # first search for _ in the input names
    match_ <- unlist(gregexpr("\\_", calaIn[i]))
    matchNum <- unlist(gregexpr("[[:digit:]]", calaIn[i]))

    if (matchNum[1] > 0) {
      temp$PARAMETER <- paste("COPE_RATIO", substring(calaIn[i], matchNum[1], nchar(calaIn[i])), sep = "_")
    } else {
      temp$PARAMETER <- paste0("COPE_RATIO", substring(calaIn[i], max(match_), nchar(calaIn[i])))
    }

    temp <- temp[, c(sampID, "PARAMETER", "RESULT")]
    metsOut <- rbind(metsOut, temp)
  }

  return(metsOut)
}
