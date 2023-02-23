#' NRSA Benthic Macroinvertebrate Taxalist
#'
#' A dataset containing the NRSA 2013-14 benthic macroinvertebrate
#' taxalist, including taxonomic and autecology information.
#'
#' @name bentTaxa_nrsa
#' @format A data frame with 1373 observations on the following 19 variables:
#'   \describe{
#'   \item{PUBLICATION_DATE}{date the dataset was published.}
#'   \item{TAXA_ID}{the taxa ID.}
#'   \item{TARGET_TAXON}{name of the taxon at the level used in NRSA.}
#'   \item{NON_TARGET}{Taxon considered non-target for NRSA but identified in
#'   sample (Y).}
#'   \item{PHYLUM}{taxonomic phylum name.}
#'   \item{CLASS}{taxonomic class name.}
#'   \item{ORDER}{taxonomic order name.}
#'   \item{FAMILY}{taxonomic family name.}
#'   \item{SUBFAMILY}{taxonomic subfamily name.}
#'   \item{TRIBE}{taxonomic tribe name.}
#'   \item{GENUS}{taxonomic genus name.}
#'   \item{FFG}{Functional feeding group trait (CF=collector-filterer,
#'   CG=collector-gatherer,PA=parasite,PI=Piercer,PR=predator,SC=scraper,SH=shredder).}
#'   \item{FFG_WSA}{Functional feeding group trait as used in WSA and NRSA metric
#'   calculations (CF=collector-filterer,CG=collector-gatherer,PA=parasite,PI=Piercer,
#'   PR=predator,SC=scraper,SH=shredder).}
#'   \item{HABIT}{Behavioral habit trait (AT=attached,BU=burrower,CB=climber,
#'   CN=clinger,DV=diver,PK=planktonic,SK=skater,SP=sprawler,SW=swimmer).}
#'   \item{HABIT_WSA}{Behavioral habit trait as used in WSA and NRSA calculations
#'   (AT=attached,BU=burrower,CB=climber,CN=clinger,DV=diver,PK=planktonic,SK=skater,
#'   SP=sprawler,SW=swimmer).}
#'   \item{PTV}{Pollution Tolerant Value.}
#'   \item{PTV_WSA}{Pollution Tolerance Values used in WSA and NRSA calculations.}
#'   \item{VOLTINISM}{Number of generations per year}
#'   \item{ITISTSN}{ITIS taxonomic serial number}

#' }
#' @note This dataset is the taxa list used for calculating metrics and indices used
#' in the NRSA 2008-9 and 2013-14 surveys. It is the default taxalist for
#' the benthic NRSA-specific functions.
#' @keywords datasets
#' @examples
#' data(bentTaxa_nrsa)
#' head(bentTaxa_nrsa)
#'
"bentTaxa_nrsa"


#' NLA Benthic Macroinvertebrate Taxalist
#'
#' A dataset containing the NLA 2012 benthic macroinvertebrate
#' taxalist, including taxonomic and autecology information.
#'
#' @name bentTaxa_nla
#' @format A data frame with 711 observations on the following 14 variables:
#'   \describe{
#'   \item{PUBLICATION_DATE}{date the dataset was published.}
#'   \item{TAXA_ID}{the taxa ID.}
#'   \item{TARGET_TAXON}{name of the taxon at the level used in NRSA.}
#'   \item{NON_TARGET}{Taxon considered non-target for NRSA but identified in
#'   sample (Y).}
#'   \item{PHYLUM}{taxonomic phylum name.}
#'   \item{CLASS}{taxonomic class name.}
#'   \item{ORDER}{taxonomic order name.}
#'   \item{FAMILY}{taxonomic family name.}
#'   \item{SUBFAMILY}{taxonomic subfamily name.}
#'   \item{TRIBE}{taxonomic tribe name.}
#'   \item{GENUS}{taxonomic genus name.}
#'   \item{FFG}{Functional feeding group trait (CF=collector-filterer,
#'   CG=collector-gatherer,PA=parasite,PI=Piercer,PR=predator,SC=scraper,SH=shredder).}
#'   \item{HABIT}{Behavioral habit trait (AT=attached,BU=burrower,CB=climber,
#'   CN=clinger,DV=diver,PK=planktonic,SK=skater,SP=sprawler,SW=swimmer).}
#'   \item{PTV}{Pollution Tolerant Value.}

#' }
#' @note This dataset is the taxa list used for calculating metrics and indices used
#' in the NLA 2007 and 2012 surveys. It is the default taxalist for
#' the benthic NRSA-specific functions.
#' @keywords datasets
#' @examples
#' data(bentTaxa_nla)
#' head(bentTaxa_nla)
#'
"bentTaxa_nla"

#' Example Benthic Invertebrate Counts
#'
#' A dataset containing benthic invertebrate count data for use in invertMet() example.
#'
#' @name bentEx
#' @format A data frame with 92 rows and the following 6 columns:
#'   \describe{
#'   \item{UID}{unique site visit ID.}
#'   \item{SAMPLE_TYPE}{sampling method used to collect sample,
#'   BERWW-wadeable reachwide sampling.}
#'   \item{SAMPLE_CAT}{identifier of (P)rimary or (D)uplicate sample for
#'   visit and SAMPLE_TYPE.}
#'   \item{TAXA_ID}{the taxa ID number corresponding to bentTaxa list.}
#'   \item{TOTAL}{number of individuals counted in sample for a given taxon.}
#'   \item{IS_DISTINCT}{indicator variable for distinctness in sample (0/1).}
#'   }
#' @note This is just a very small subset of benthic count data for NRSA 2008-2009
#'    for example purposes only.
#' @examples
#'    data(bentEx)
#'    head(bentEx)
#' @keywords datasets
"bentEx"

#' Example Distinctness Input Benthic Dataset
#'
#' A dataset containing benthic invertebrate count data for use in assignDistinct() example.
#'
#' @name distEx
#' @format A data frame with 92 observations on the following 11 variables:
#'   \describe{
#'   \item{TAXA_ID}{the taxa ID number corresponding to fishTaxa list.}
#'   \item{UID}{unique site visit ID.}
#'   \item{SAMPLE_TYPE}{sampling method used to collect sample, BERWW-wadeable
#'   reachwide sampling.}
#'   \item{TARGET_TAXON}{name of the taxon at the level used in NRSA.}
#'   \item{TOTAL}{number of individuals counted in sample for a given taxon.}
#'   \item{PHYLUM}{taxonomic phylum name.}
#'   \item{CLASS}{taxonomic class name.}
#'   \item{ORDER}{taxonomic order name.}
#'   \item{FAMILY}{taxonomic family name.}
#'   \item{GENUS}{taxonomic genus name.}
#'   \item{SAMPID}{variable to uniquely identify samples, here the
#'   concatenation of UID and SAMPLE_TYPE}
#'   }
#' @note This is just a very small subset of benthic count data for NRSA 2008-2009
#'   for example purposes only.
#' @examples
#'   data(distEx)
#'   head(distEx)
#' @keywords datasets
"distEx"

#' Example Fish Count Dataset
#'
#' A dataset containing fish count data for use in fishMet() example.
#'
#' @name fishEx
#' @format A data frame with 22 observations on the following 6 variables.
#'   \describe{
#'   \item{UID}{unique site visit ID.}
#'   \item{TAXA_ID}{the taxa ID number corresponding to fishTaxa list.}
#'   \item{FINAL_CT}{number of fish captured for a given taxon.}
#'   \item{ANOM_CT}{number of anomalies observed for a given taxon.}
#'   \item{IS_DISTINCT}{indicator variable for distinctness in sample (0/1).}
#'   \item{NON_NATIVE}{Indicator variable of site-specific non-native status (Y/N).}
#'   }
#' @note This is just a very small subset of fish count data for NRSA 2008-2009
#'   for example purposes only.
#' @examples
#'   data(fishEx)
#'   head(fishEx)
#' @keywords datasets
"fishEx"

#' Fish Taxa
#'
#' A dataset containing fish taxa.
#'
#' @name fishTaxa
#' @format A data frame with 851 rows and the following 16 columns:
#'   \describe{
#'   \item{PUBLICATION_DATE}{date the dataset was published.}
#'   \item{TAXA_ID}{the taxa ID.}
#'   \item{FAM_OR_CLS}{taxonomic family or class name.}
#'   \item{FAMILY}{family name.}
#'   \item{FINAL_NAME}{species common name.}
#'   \item{GENUS}{taxonomic genus name.}
#'   \item{HABITAT_NRSA}{HABITAT preference of taxon used in NRSA calculations.}
#'   \item{HERP}{Taxon is either reptile (R) or amphibian (A).}
#'   \item{ITISTSN}{ITIS taxonomic serial number.}
#'   \item{MIGR_NRSA}{Migratory trait of taxon (Y/N).}
#'   \item{REPROD_NRSA}{Reproductive strategy of taxon (C=clean,coarse (lithophil),
#'   D=drifter,G=guarder,O=other).}
#'   \item{TEMP_NRSA}{Temperature preference trait used in NRSA (WM=warm, CL=cool, CD=cold).}
#'   \item{TOL_VAL_EMAPW}{Species tolerance values based on EMAP-West data.}
#'   \item{TOLERANCE_NRSA}{Tolerance categories used in NRSA based on WEMAP, USGS,
#'   Region 7 (S=sensitive/intolerant,I=intermediate,T=tolerant).}
#'   \item{TROPHIC_NRSA}{Trophic guild as used in NRSA (C=carnivore,I=invertivore,
#'   H=herbivore,O=omnivore).}
#'   \item{VEL_NRSA}{Velocity preference as used in NRSA (R=rheophil,P=pool,O=other).}
#'   }
#' @note Unless a data frame name is provided for argument inTaxa to the fishMet
#'   function, this dataset provides the inTaxa data frame.
#' @examples
#'   data(fishTaxa)
#'   head(fishTaxa)
#' @keywords datasets
"fishTaxa"

#' Zooplankton Taxa
#'
#' A dataset containing zooplankton taxa.
#'
#' @name zoopTaxa
#' @format A data frame with 581 rows and the following 18 columns:
#' \describe{
#' \item{PUBLICATION_DATE}{date the dataset was published.}
#' \item{TAXA_ID}{the taxa ID.}
#' \item{TARGET_TAXON}{name of the taxon at the level used in NRSA.}
#' \item{PHYLUM}{taxonomic phylum name.}
#' \item{CLASS}{taxonomic class name.}
#' \item{SUBCLASS}{taxonomic subclass name.}
#' \item{ORDER}{taxonomic order name.}
#' \item{SUBORDER}{taxonomic suborder name.}
#' \item{FAMILY}{taxonomic family name.}
#' \item{SUBFAMILY}{taxonomic subfamily name.}
#' \item{GENUS}{taxonomic genus name.}
#' \item{SPECIES}{taxonomic species name.}
#' \item{SUBSPECIES}{taxonomic subspecies name.}
#' \item{FFG}{Functional feeding group trait (OMNI = Omnivore,
#' HERB = Herbivore, PARA = Parasite, FILT = Filterer,
#' PRED = Predator, UNK = Unknown)}
#' \item{NON_TARGET}{Non-target taxon for NLA purposes (Y or NA)}
#' \item{NON_NATIVE}{Non-native status. If value is ALL, non-native
#' in all states; if value is a list of states, states where taxon is
#' non-native; if missing, taxon considered native.}
#' \item{NET_SIZECLS_NEW}{Mesh size class assigned to taxon (COARSE/FINE)
#' for NLA}
#' \item{CLADOCERA_SIZE}{Size class of cladoceran taxon (LARGE/SMALL)}
#' }
#' @note Unless a data frame name is provided for argument inTaxa to the
#' zooplankton metric functions, this dataset provides the inTaxa data frame.
#' @examples
#'   data(zoopTaxa)
#'   head(zoopTaxa)
#' @keywords datasets
"zoopTaxa"


