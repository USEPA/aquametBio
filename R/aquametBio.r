#' aquametBio: Calculate fish and benthic macroinvertebrate metrics used in NRSA and NLA
#'
#' @docType package
#' @name aquametBio
#'
#' @importFrom Hmisc "%nin%" capitalize
#' @importFrom gtools smartbind
#' @importFrom plyr ddply mutate summarize summarise rename revalue
#' @importFrom dplyr filter select "%>%" group_by
#' @importFrom reshape2 dcast melt
#' @importFrom stringr str_detect
#' @importFrom stats aggregate approx
#' @importFrom utils head
#' @keywords package
#' @title aquametBio
#'
#' @keywords package
#' @title aquametBio

if(getRversion() >= "3.0") utils::globalVariables(c('IS_DISTINCT','FINAL_CT'
,'domN','TOTSUM','TOTSUM','SAMP_ID','WSAREA','TOTLNIND','MMI_FISH'
,'gf','fp','FISH_MMI_COND','NON_TARGET','fishTaxa','SAMPID','TAXA_ID'
,'DOM3PIND','DOM1PIND','DOM5PIND','FAMILY','CHIRDOM1PIND','CHIRDOM3PIND'
,'CHIRDOM5PIND','value','FFG','ORDER','TOTLNTAX','HABIT','CLASS'
,'PHYLUM','SUBFAMILY','TRIBE','ORTHNIND','CHIRNIND','ORTHCHIRPIND'
,'PTV','TRAIT','variable','ANOM_CT','ANOMPIND','PARAMETER','Freq'
,'RESULT','int','slope','LWSAREA','RESULT_WS','NONNATIVE','ALIEN'
,'NAT','VELOCITY','MIGRATORY','REPROD','TEMP','NAT_TOTLNIND'
,'NAT_TOTLNTAX','GENUS','NAME','TOLERANCE','INTL','NTOL','HABITAT'
,'TROPHIC','TOL_VAL','bentTaxa_nla','ECO_BIO','MMI_BENT'
,'bentTaxa_nrsa','ECO9','bentTaxa','TARGET_TAXON','TOTAL'
,'SUMCT'))




