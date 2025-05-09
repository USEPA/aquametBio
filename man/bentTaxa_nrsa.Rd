% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.r
\docType{data}
\name{bentTaxa_nrsa}
\alias{bentTaxa_nrsa}
\title{NRSA Benthic Macroinvertebrate Taxalist}
\format{
A data frame with 1373 observations on the following 19 variables:
  \describe{
  \item{PUBLICATION_DATE}{date the dataset was published.}
  \item{TAXA_ID}{the taxa ID.}
  \item{TARGET_TAXON}{name of the taxon at the level used in NRSA.}
  \item{NON_TARGET}{Taxon considered non-target for NRSA but identified in
  sample (Y).}
  \item{PHYLUM}{taxonomic phylum name.}
  \item{CLASS}{taxonomic class name.}
  \item{ORDER}{taxonomic order name.}
  \item{FAMILY}{taxonomic family name.}
  \item{SUBFAMILY}{taxonomic subfamily name.}
  \item{TRIBE}{taxonomic tribe name.}
  \item{GENUS}{taxonomic genus name.}
  \item{FFG}{Functional feeding group trait (CF=collector-filterer,
  CG=collector-gatherer,PA=parasite,PI=Piercer,PR=predator,SC=scraper,SH=shredder).}
  \item{FFG_WSA}{Functional feeding group trait as used in WSA and NRSA metric
  calculations (CF=collector-filterer,CG=collector-gatherer,PA=parasite,PI=Piercer,
  PR=predator,SC=scraper,SH=shredder).}
  \item{HABIT}{Behavioral habit trait (AT=attached,BU=burrower,CB=climber,
  CN=clinger,DV=diver,PK=planktonic,SK=skater,SP=sprawler,SW=swimmer).}
  \item{HABIT_WSA}{Behavioral habit trait as used in WSA and NRSA calculations
  (AT=attached,BU=burrower,CB=climber,CN=clinger,DV=diver,PK=planktonic,SK=skater,
  SP=sprawler,SW=swimmer).}
  \item{PTV}{Pollution Tolerant Value.}
  \item{PTV_WSA}{Pollution Tolerance Values used in WSA and NRSA calculations.}
  \item{VOLTINISM}{Number of generations per year}
  \item{ITISTSN}{ITIS taxonomic serial number}
}
}
\usage{
bentTaxa_nrsa
}
\description{
A dataset containing the NRSA 2013-14 benthic macroinvertebrate
taxalist, including taxonomic and autecology information.
}
\note{
This dataset is the taxa list used for calculating metrics and indices used
in the NRSA 2008-9 and 2013-14 surveys. It is the default taxalist for
the benthic NRSA-specific functions.
}
\examples{
data(bentTaxa_nrsa)
head(bentTaxa_nrsa)

}
\keyword{datasets}
