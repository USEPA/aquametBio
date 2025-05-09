% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcZoopCopeMetrics.r
\name{calcZoopCopeMetrics}
\alias{calcZoopCopeMetrics}
\title{Calculate zooplankton copepod ratio metrics}
\usage{
calcZoopCopeMetrics(indata, sampID, calaIn, cyclIn, cladIn)
}
\arguments{
\item{indata}{Input data frame, containing \emph{sampID}
as variable identifying unique samples, as well as the
metrics listed in \emph{calaIn}, \emph{cyclIn}, and
\emph{cladIn}.}

\item{sampID}{A character vector containing the names of all
variables in \emph{indata} that specify a unique sample.}

\item{calaIn}{A vector of the names of input Calanoida metrics, in
the same order as for cladIn and cyclIn metrics.}

\item{cyclIn}{A vector of the names of input Cyclopoida metrics, in
the same order as for calaIn and cladIn metrics.}

\item{cladIn}{A vector of the names of input Cladocera metrics, in
the same order as for calaIn and cyclIn metrics.}
}
\value{
A data frame with \emph{sampID} variables and the metric
containing the copepod ratio metrics

\href{https://github.com/USEPA/aquametBio/blob/main/inst/NLA_Zooplankton_Metric_Descriptions.pdf}{NLA_Zooplankton_Metric_Descriptions.pdf}
}
\description{
This function calculates ratios of calanoid
copepods to the sum of cyclopoid copepods and cladocerans
}
\author{
Karen Blocksom \email{Blocksom.Karen@epa.gov}
}
