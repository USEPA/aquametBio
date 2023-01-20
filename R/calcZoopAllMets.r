#' @export
#' @title Calculate all zooplankton metrics
#'
#' @description This function calculates all zooplankton metrics
#' for a given input data frame using the inputs provided. These
#' inputs include the names of count, biomass, and density variables,
#' as well as a variable indicating distinctness of each taxon.
#' The default value for each argument is that used for NLA.
#'
#' @param indata A data frame containing, at minimum, the variables
#' specified in the arguments for sampID, is_distinct, ct, biomass,
#' density, ct_sub, biomass_sub, taxa_id, nonnative. The dataset
#' should contain only data from the combined sample type (ZONW for NLA).
#' @param inCoarse A data frame for coarse mesh zooplankton sample data,
#' containing the variables
#' specified in the arguments for sampID, is_distinct, ct, biomass,
#' density, ct_sub, biomass_sub, taxa_id, nonnative. The dataset
#' should contain only data from the combined sample type (ZOCN for NLA)
#' @param inFine A data frame for fine mesh zooplankton sample data,
#' containing the variables
#' specified in the arguments for sampID, is_distinct, ct, biomass,
#' density, ct_sub, biomass_sub, taxa_id, nonnative. The dataset
#' should contain only data from the combined sample type (ZOFN for NLA)
#' @param inTaxa a data frame containing taxonomic information,
#' including the variables PHYLUM, CLASS, SUBCLASS, ORDER, SUBORDER,
#' FAMILY, following the same taxonomy as the default taxa list,
#' \emph{zoopTaxa}. There should also be autecology traits with names
#' that match those for the arguments \emph{ffg}, \emph{clad_size},
#' \emph{net_size}, \emph{genus}, and \emph{family}.
#' @param sampID A character vector containing the names of all
#' variables in indf that specify a unique sample. If not specified,
#' the default is \emph{UID}
#' @param is_distinct A string with the name of the distinctness variable,
#' which is assumed to have only values of 0 or 1. If not specified,
#' the default is \emph{IS_DISTINCT}.
#' @param ct A string with the name of the count variable. This
#' variable should be numeric or able to be converted to numeric.
#' @param biomass A string with the name of the biomass variable. This
#' variable should be numeric or able to be converted to numeric.
#' @param density A string with the name of the density variable. This
#' variable should be numeric or able to be converted to numeric.
#' @param is_distinct_sub A string with the name of the distinctness
#' variable for the combined subsample,
#' which is assumed to have only values of 0 or 1. If not specified,
#' the default is \emph{IS_DISTINCT_300}.
#' @param ct_sub A string with the name of the count variable for
#' the combined subsample.
#' @param biomass_sub A string with the name of the biomass
#' variable for the combined subsample.
#' @param sub_mod A string with the modifier used in the names
#' for the \emph{ct_sub} and \emph{biomass_sub} variables. The
#' default is '300' because this is the subsample size used for
#' individual coarse and fine mesh samples for NLA.
#' @param taxa_id A string with the name of the taxon ID variable
#' in \emph{inCts} that matches that in \emph{inTaxa}. The default
#' value is \emph{TAXA_ID}.
#' @param nonnative A string with the name of the numeric variable
#' indicating a non-native taxon. A value of 1 indicates that a
#' taxon is non-native, and 0 native. All other values will be ignored.
#' @param genus A string with the name of the variable containing
#' genus name in the \emph{intaxa} data frame. Blank or missing
#' genus names will be dropped.
#' @param family A string with the name of the variable containing
#' family name in the \emph{intaxa} data frame. Blank or missing
#' family names will be dropped.
#' @param ffg A string with the name of the functional feeding group
#' variable in \emph{inTaxa}. Values used
#' in calculations include FILT, HERB, OMNI, PARA, PRED, representing
#' filterer, herbivore, omnivore, parasite, and predator,
#' respectively.
#' @param clad_size A string with the name of the variable in
#' \emph{inTaxa} indicating the size class of the cladoceran
#' taxon. Valid values are LARGE and SMALL.
#' @param net_size A string with the name of the variable in
#' \emph{inTaxa} indicating the net size class of a taxon.
#' Valid values are COARSE and FINE.
calcZoopAllMets.r <- function(indata, inCoarse, inFine,
                              inTaxa, sampID, is_distinct,
                              ct = 'COUNT', biomass = 'BIOMASS',
                              density = 'DENSITY',
                              is_distinct_sub = 'IS_DISTINCT_300',
                              ct_sub = 'COUNT_300',
                              biomass_sub = 'BIOMASS_300',
                              sub_mod = '300', taxa_id = 'TAXA_ID',
                              nonnative = 'NON_NATIVE', genus = 'GENUS',
                              family = 'FAMILY', ffg = 'FFG',
                              clad_size = 'CLADOCERA_SIZE',
                              net_size = 'NET_SIZECLS_NEW'){

  indata <- as.data.frame(indata)

  necVars <- c(sampID, is_distinct, ct, biomass, density, ct_sub,
               biomass_sub, taxa_id, nonnative)
  if(any(necVars %nin% names(indata))){
    msgTraits <- which(necVars %nin% names(indata))
    print(paste("Missing variables in input data frame:",
                paste(necVars[msgTraits], collapse=',')))
    return(NULL)
  }else{
    indata[, c(ct, biomass, ct_sub, density, biomass_sub, is_distinct)] <- lapply(indata[, c(ct, biomass, ct_sub, density, biomass_sub, is_distinct)], as.numeric)
  }

  necTaxVars <- c(taxa_id, ffg, clad_size, net_size,
                  'PHYLUM', 'CLASS', 'SUBCLASS', 'ORDER',
                  'SUBORDER', 'FAMILY', genus, family)
  if(any(necTaxVars %nin% names(inTaxa))){
    msgTraits <- which(necTaxVars %nin% names(inTaxa))
    print(paste("Missing variables in input taxalist:",
                paste(necTaxVars[msgTraits], collapse=',')))
    return(NULL)
  }

  baseMets.full <- calcZoopBaseMetrics(indata, sampID, is_distinct,
                                       ct, biomass, density,
                                       inTaxa, taxa_id,
                                       ffg, clad_size,
                                       net_size,
                                       nativeMetrics = FALSE)

  indata.nat <- subset(indata, eval(as.name(nonnative))==0)

  baseMets.nat <- calcZoopBaseMetrics(indata.nat, sampID, is_distinct,
                                       ct, biomass, density,
                                       inTaxa, taxa_id,
                                       ffg, clad_size,
                                       net_size,
                                       nativeMetrics = TRUE)

  baseMets.300 <- calcZoopBaseMetrics(indata.nat, sampID, is_distinct_sub,
                                      ct_sub, biomass_sub,
                                      inTaxa, taxa_id,
                                      ffg, clad_size,
                                      net_size,
                                      nativeMetrics = FALSE)

  baseMets.300.nat <- calcZoopBaseMetrics(indata.nat, sampID, is_distinct_sub,
                                          ct_sub, biomass_sub,
                                          inTaxa, taxa_id,
                                          ffg, clad_size,
                                          net_size,
                                          nativeMetrics = TRUE)


  # Copepod ratio metrics - these require metrics calculated above
  copeMetsIn <- rbind(subset(baseMets.full, PARAMETER %in% c('CALAN_NIND',
                                                    'CALAN_BIO', 'CALAN_DEN',
                                                    'CYCLOP_NIND', 'CYCLOP_BIO',
                                                    'CYCLOP_DEN', 'CLAD_NIND',
                                                    'CLAD_BIO', 'CLAD_DEN')),
                      subset(baseMets.300, PARAMETER %in% c('CALAN300_NIND',
                                                    'CALAN300_BIO',
                                                    'CYCLOP300_NIND', 'CYCLOP300_BIO',
                                                    'CLAD300_NIND', 'CLAD300_BIO')))


  copeRatMets <- calcZoopCopeMetrics(copeMetsIn, sampID,
                                     c('CALAN_NIND', 'CALAN_BIO', 'CALAN_DEN',
                                       'CALAN300_NIND', 'CALAN300_BIO'),
                                     c('CYCLOP_NIND', 'CYCLOP_BIO', 'CYCLOP_DEN',
                                       'CYCLOP300_NIND', 'CYCLOP300_BIO'),
                                     c('CLAD_NIND', 'CLAD_BIO', 'CLAD_DEN',
                                       'CLAD300_NIND', 'CLAD300_BIO'))

# Diversity metrics
  divMets.full <- calcZoopDivMetrics(indata, sampID, is_distinct,
                                     ct, biomass, density)

  divMets.300 <- calcZoopDivMetrics(indata, sampID, is_distinct_sub,
                                    ct_sub, biomass_sub)

  # Now merge taxa with indata to subset to rotifers
  indata.taxa <- merge(indata, inTaxa, by = taxa_id) |>
    subset(select = c(sampID, ct, ct_sub, taxa_id, is_distinct, is_distinct_sub,
                      'SUBORDER', 'SUBCLASS', 'PHYLUM'))

  divMets.full.rot <- calcZoopDivMetrics(subset(indata.taxa, PHYLUM=='ROTIFERA'),
                                         sampID, is_distinct, ct)

  divMets.full.clad <- calcZoopDivMetrics(subset(indata.taxa, SUBORDER=='CLADOCERA'),
                                         sampID, is_distinct, ct)

  divMets.full.cope <- calcZoopDivMetrics(subset(indata.taxa, SUBCLASS=='COPEPODA'),
                                         sampID, is_distinct, ct)

  divMets.300.rot <- calcZoopDivMetrics(subset(indata.taxa, PHYLUM=='ROTIFERA'),
                                         sampID, is_distinct_sub, ct_sub)

  divMets.300.clad <- calcZoopDivMetrics(subset(indata.taxa, SUBORDER=='CLADOCERA'),
                                          sampID, is_distinct_sub, ct_sub)

  divMets.300.cope <- calcZoopDivMetrics(subset(indata.taxa, SUBCLASS=='COPEPODA'),
                                          sampID, is_distinct_sub, ct_sub)
# Dominance metrics
  dom.full <- calcZoopDomMetrics(indata, sampID, is_distinct,
                                  valsIn = c(ct, biomass, density),
                                  valsOut = c('PIND', 'PBIO', 'PDEN'),
                                  taxa_id,
                                  subgrp = NULL)

  dom.300 <- calcZoopDomMetrics(indata, sampID, is_distinct_sub,
                                 valsIn = c(ct_sub, biomass_sub),
                                 valsOut = c('300_PIND', '300_PBIO', '300_PDEN'),
                                 taxa_id,
                                 subgrp = NULL)

  dom.full.rot <- calcZoopDomMetrics(subset(indata.taxa, PHYLUM=='ROTIFERA'),
                                     sampID, is_distinct,
                                     valsIn = c(ct, biomass, density),
                                     valsOut = c('PIND', 'PBIO', 'PDEN'),
                                     taxa_id, subgrp = 'ROT')

  dom.full.clad <- calcZoopDomMetrics(subset(indata.taxa, SUBORDER=='CLADOCERA'),
                                     sampID, is_distinct,
                                     valsIn = c(ct, biomass, density),
                                     valsOut = c('PIND', 'PBIO', 'PDEN'),
                                     taxa_id, subgrp = 'CLAD')

  dom.full.cope <- calcZoopDomMetrics(subset(indata.taxa, SUBCLASS=='COPEPODA'),
                                     sampID, is_distinct,
                                     valsIn = c(ct, biomass, density),
                                     valsOut = c('PIND', 'PBIO', 'PDEN'),
                                     taxa_id, subgrp = 'COPE')

  dom.300.rot <- calcZoopDomMetrics(subset(indata.taxa, PHYLUM=='ROTIFERA'),
                                     sampID, is_distinct_sub,
                                     valsIn = c(ct_sub, biomass_sub),
                                     valsOut = c('PIND', 'PBIO'),
                                     taxa_id, subgrp = 'ROT')
  dom.300.rot$PARAMETER <- gsub('PIND', '300_PIND', dom.300.rot$PARAMETER)
  dom.300.rot$PARAMETER <- gsub('PBIO', '300_PBIO', dom.300.rot$PARAMETER)

  dom.300.clad <- calcZoopDomMetrics(subset(indata.taxa, SUBORDER=='CLADOCERA'),
                                     sampID, is_distinct_sub,
                                     valsIn = c(ct_sub, biomass_sub),
                                      valsOut = c('PIND', 'PBIO'),
                                      taxa_id, subgrp = 'CLAD')
  dom.300.clad$PARAMETER <- gsub('PIND', '300_PIND', dom.300.clad$PARAMETER)
  dom.300.clad$PARAMETER <- gsub('PBIO', '300_PBIO', dom.300.clad$PARAMETER)

  dom.300.cope <- calcZoopDomMetrics(subset(indata.taxa, SUBCLASS=='COPEPODA'),
                                     sampID, is_distinct_sub,
                                     valsIn = c(ct_sub, biomass_sub),
                                     valsOut = c('PIND', 'PBIO'),
                                     taxa_id, subgrp = 'COPE')
  dom.300.cope$PARAMETER <- gsub('PIND', '300_PIND', dom.300.cope$PARAMETER)
  dom.300.cope$PARAMETER <- gsub('PBIO', '300_PBIO', dom.300.cope$PARAMETER)

  # Richness metrics
  richMets.full <- calcZoopRichnessMetrics(indata, sampID,
                                           distVars=c(is_distinct, is_distinct_sub),
                                           nonnative, inTaxa,
                                           taxa_id, genus,
                                           family, prefix = c('', '300'))

  # Totals - these are all in wide format so combine and then melt
  totMets.full <- calcZoopTotals(indata, sampID, is_distinct,
                                 c(ct, biomass, density),
                                 c('TOTL_NIND', 'TOTL_BIO', 'TOTL_DEN'),
                                 'TOTL_NTAX')

  totMets.nat <- calcZoopTotals(indata.nat, sampID, is_distinct,
                                c(ct, biomass, density),
                                c('TOTL_NAT_NIND', 'TOTL_NAT_BIO', 'TOTL_NAT_DEN'),
                                'TOTL_NAT_NTAX')

  totMets.300 <- calcZoopTotals(indata, sampID, is_distinct_sub,
                                 c(ct_sub, biomass_sub),
                                 c('TOTL300_NIND', 'TOTL300_BIO'),
                                 'TOTL300_NTAX')

  totMets.300.nat <- calcZoopTotals(indata.nat, sampID, is_distinct_sub,
                                c(ct_sub, biomass_sub),
                                c('TOTL300_NAT_NIND', 'TOTL300_NAT_BIO'),
                                'TOTL300_NAT_NTAX')

  # Also calculate totals by coarse and fine sample types
  totMets.zocn.full <- calcZoopTotals(inCoarse, sampID, is_distinct,
                                      c(ct, biomass, density),
                                      c('ZOCN_NIND', 'ZOCN_BIO', 'ZOCN_DEN'),
                                      'ZOCN_NTAX')

  totMets.zocn.nat <- calcZoopTotals(subset(inCoarse, eval(as.name(nonnative))==0),
                                     sampID, is_distinct,
                                     c(ct, biomass, density),
                                     c('ZOCN_NAT_NIND', 'ZOCN_NAT_BIO', 'ZOCN_NAT_DEN'),
                                     'ZOCN_NAT_NTAX')

  totMets.zocn.300 <- calcZoopTotals(inCoarse, sampID, is_distinct_sub,
                                     c(ct_sub, biomass_sub),
                                     c('ZOCN300_NIND', 'ZOCN300_BIO'),
                                     'ZOCN300_NTAX')

  totMets.zocn.300.nat <- calcZoopTotals(subset(inCoarse, eval(as.name(nonnative))==0),
                                         sampID, is_distinct_sub,
                                         c(ct_sub, biomass_sub),
                                         c('ZOCN300_NAT_NIND', 'ZOCN300_NAT_BIO'),
                                         'ZOCN300_NAT_NTAX')

  totMets.zofn.full <- calcZoopTotals(inFine, sampID, is_distinct,
                                      c(ct, biomass, density),
                                      c('ZOFN_NIND', 'ZOFN_BIO', 'ZOFN_DEN'),
                                      'ZOFN_NTAX')

  totMets.zofn.nat <- calcZoopTotals(subset(inFine, eval(as.name(nonnative))==0),
                                     sampID, is_distinct,
                                     c(ct, biomass, density),
                                     c('ZOFN_NAT_NIND', 'ZOFN_NAT_BIO', 'ZOFN_NAT_DEN'),
                                     'ZOFN_NAT_NTAX')

  totMets.zofn.300 <- calcZoopTotals(inFine, sampID, is_distinct_sub,
                                     c(ct_sub, biomass_sub),
                                     c('ZOFN300_NIND', 'ZOFN300_BIO'),
                                     'ZOFN300_NTAX')

  totMets.zofn.300.nat <- calcZoopTotals(subset(inFine, eval(as.name(nonnative))==0),
                                         sampID, is_distinct_sub,
                                         c(ct_sub, biomass_sub),
                                         c('ZOFN300_NAT_NIND', 'ZOFN300_NAT_BIO'),
                                         'ZOFN300_NAT_NTAX')

  zonwMets.all <- merge(totMets.full, totMets.nat, by = sampID, all.x=TRUE) |>
    merge(totMets.300, by = sampID, all.x=TRUE) |>
    merge(totMets.300.nat, by = sampID, all.x=TRUE)

  zocnTot.all <- merge(totMets.zocn.full, totMets.zocn.nat, by = sampID, all.x=TRUE) |>
    merge(totMets.zocn.300, by = sampID, all.x=TRUE) |>
    merge(totMets.zocn.300.nat, by = sampID, all.x=TRUE)

  zofnTot.all <- merge(totMets.zofn.full, totMets.zofn.nat, by = sampID, all.x=TRUE) |>
    merge(totMets.zofn.300, by = sampID, all.x=TRUE) |>
    merge(totMets.zofn.300.nat, by = sampID, all.x=TRUE)

  totMets.all <- merge(zonwMets.all, zocnTot.all, by = sampID, all.x=TRUE) |>
    merge(zofnTot.all, by = sampID, all.x=TRUE)

  totMets.all.long <- reshape(totMets.all, idvar = sampID, direction = 'long',
            varying = names(totMets.all)[!(names(totMets.all) %in% sampID)],
            timevar = 'PARAMETER', v.names = 'RESULT',
            times = names(totMets.all)[!(names(totMets.all) %in% sampID)])

  # Native metrics
  natMets.full <- calcZoopNativeMetrics(matchSums, c('UID', 'SAMPLE_TYPE'),
                        inputNative = c('TOTL_NAT_NTAX', 'TOTL300_NAT_NTAX',
                                        'ZOCN_NAT_NTAX', 'ZOCN300_NAT_NTAX',
                                        'ZOFN_NAT_NTAX', 'ZOFN300_NAT_NTAX',
                                        'TOTL_NAT_DEN', 'TOTL_NAT_BIO',
                                        'TOTL_NAT_NIND', 'TOTL300_NAT_NIND',
                                        'TOTL300_NAT_BIO', 'ZOCN_NAT_NIND',
                                        'ZOCN300_NAT_NIND', 'ZOCN_NAT_BIO',
                                        'ZOCN300_NAT_BIO', 'ZOCN_NAT_DEN',
                                        'ZOFN_NAT_NIND', 'ZOFN300_NAT_NIND',
                                        'ZOFN_NAT_BIO', 'ZOFN300_NAT_BIO'),
                        inputTotals = c('TOTL_NTAX', 'TOTL300_NTAX',
                                        'ZOCN_NTAX', 'ZOCN300_NTAX',
                                        'ZOFN_NTAX', 'ZOFN300_NTAX',
                                        'TOTL_DEN',
                                        'TOTL_BIO', 'TOTL_NIND',
                                        'TOTL300_NIND', 'TOTL300_BIO',
                                        'ZOCN_NIND', 'ZOCN300_NIND',
                                        'ZOCN_BIO', 'ZOCN300_BIO',
                                        'ZOCN_DEN', 'ZOFN_NIND',
                                        'ZOFN300_NIND',
                                        'ZOFN_BIO', 'ZOFN300_BIO'))

  # Now combine all metrics together
  allMets <- rbind(baseMets.full, baseMets.nat, baseMets.300,
                   baseMets.300.nat, copeRatMets,
                   divMets.full, divMets.300,
                   divMets.full.rot, divMets.full.clad,
                   divMets.full.cope, divMets.300.rot,
                   divMets.300.clad, divMets.300.cope,
                   dom.full, dom.300, dom.full.rot,
                   dom.full.clad, dom.full.cope,
                   dom.300.rot, dom.300.clad,
                   dom.300.cope, richMets.full,
                   totMets.all.long, natMets.full)
}
