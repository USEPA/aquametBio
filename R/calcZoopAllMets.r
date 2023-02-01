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
calcZoopAllMets <- function(indata, inCoarse, inFine,
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
    indata[, c(ct, biomass, ct_sub, density, biomass_sub, is_distinct, is_distinct_sub)] <- lapply(indata[, c(ct, biomass, ct_sub, density, biomass_sub, is_distinct, is_distinct_sub)], as.numeric)
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

  baseMets.sub <- calcZoopBaseMetrics(indata, sampID, is_distinct_sub,
                                      ct_sub, biomass_sub, density = NULL,
                                      inTaxa, taxa_id,
                                      ffg, clad_size,
                                      net_size,
                                      nativeMetrics = FALSE)

  baseMets.sub$PARAMETER <- gsub("\\_NIND", paste0(sub_mod, "\\_NIND"), baseMets.sub$PARAMETER)
  baseMets.sub$PARAMETER <- gsub("\\_PIND", paste0(sub_mod, "\\_PIND"), baseMets.sub$PARAMETER)
  baseMets.sub$PARAMETER <- gsub("\\_NTAX", paste0(sub_mod, "\\_NTAX"), baseMets.sub$PARAMETER)
  baseMets.sub$PARAMETER <- gsub("\\_PTAX", paste0(sub_mod, "\\_PNTAX"), baseMets.sub$PARAMETER)
  baseMets.sub$PARAMETER <- gsub("\\_BIO", paste0(sub_mod, "\\_BIO"), baseMets.sub$PARAMETER)
  baseMets.sub$PARAMETER <- gsub("\\_PBIO", paste0(sub_mod, "\\_PBIO"), baseMets.sub$PARAMETER)

  baseMets.sub.nat <- calcZoopBaseMetrics(indata.nat, sampID, is_distinct_sub,
                                          ct_sub, biomass_sub, density = NULL,
                                          inTaxa, taxa_id,
                                          ffg, clad_size,
                                          net_size,
                                          nativeMetrics = TRUE)

  baseMets.sub.nat$PARAMETER <- gsub("\\_NAT\\_NIND", paste0(sub_mod, "\\_NAT\\_NIND"), baseMets.sub.nat$PARAMETER)
  baseMets.sub.nat$PARAMETER <- gsub("\\_NAT\\_PIND", paste0(sub_mod, "\\_NAT\\_PIND"), baseMets.sub.nat$PARAMETER)
  baseMets.sub.nat$PARAMETER <- gsub("\\_NAT\\_NTAX", paste0(sub_mod, "\\_NAT\\_NTAX"), baseMets.sub.nat$PARAMETER)
  baseMets.sub.nat$PARAMETER <- gsub("\\_NAT\\_PTAX", paste0(sub_mod, "\\_NAT\\_PNTAX"), baseMets.sub.nat$PARAMETER)
  baseMets.sub.nat$PARAMETER <- gsub("\\_NAT\\_BIO", paste0(sub_mod, "\\_NAT\\_BIO"), baseMets.sub.nat$PARAMETER)
  baseMets.sub.nat$PARAMETER <- gsub("\\_NAT\\_PBIO", paste0(sub_mod, "\\_NAT\\_PBIO"), baseMets.sub.nat$PARAMETER)

  print("Finished with base metrics")

  # Copepod ratio metrics - these require metrics calculated above
  copeMetsIn <- rbind(subset(baseMets.full, PARAMETER %in% c('CALAN_NIND',
                                                    'CALAN_BIO', 'CALAN_DEN',
                                                    'CYCLOP_NIND', 'CYCLOP_BIO',
                                                    'CYCLOP_DEN', 'CLAD_NIND',
                                                    'CLAD_BIO', 'CLAD_DEN')),
                      subset(baseMets.sub, PARAMETER %in% c(paste0('CALAN', sub_mod, '_NIND'),
                                                    paste0('CALAN', sub_mod, '_BIO'),
                                                    paste0('CYCLOP', sub_mod, '_NIND'),
                                                    paste0('CYCLOP', sub_mod, '_BIO'),
                                                    paste0('CLAD', sub_mod, '_NIND'),
                                                    paste0('CLAD', sub_mod, '_BIO'))))

  copeMetsIn.wide <- reshape(copeMetsIn, idvar = c(sampID), direction = 'wide',
                             timevar = 'PARAMETER', v.names = 'RESULT')

  names(copeMetsIn.wide) <- gsub("RESULT\\.", "", names(copeMetsIn.wide))


  copeRatMets <- calcZoopCopeMetrics(copeMetsIn.wide, sampID,
                                     c('CALAN_NIND', 'CALAN_BIO', 'CALAN_DEN',
                                       paste0('CALAN', sub_mod, '_NIND'),
                                       paste0('CALAN', sub_mod, '_BIO')),
                                     c('CYCLOP_NIND', 'CYCLOP_BIO', 'CYCLOP_DEN',
                                       paste0('CYCLOP', sub_mod, '_NIND'),
                                       paste0('CYCLOP', sub_mod, '_BIO')),
                                     c('CLAD_NIND', 'CLAD_BIO', 'CLAD_DEN',
                                       paste0('CLAD', sub_mod, '_NIND'),
                                       paste0('CLAD', sub_mod, '_BIO')))

  print("Finished copepod ratio merics")

# Diversity metrics
  divMets.full <- calcZoopDivMetrics(indata, sampID, is_distinct,
                                     ct, biomass, density)

  divMets.sub <- calcZoopDivMetrics(indata, sampID, is_distinct_sub,
                                    ct_sub, biomass_sub)

  divMets.sub$PARAMETER <- gsub("\\_NIND", paste0(sub_mod, "\\_NIND"), divMets.sub$PARAMETER)
  divMets.sub$PARAMETER <- gsub("\\_BIO", paste0(sub_mod, "\\_BIO"), divMets.sub$PARAMETER)

  print("Finished overall diversity metrics")

  # Now merge taxa with indata to subset to rotifers
  indata.taxa <- merge(indata, inTaxa, by = taxa_id) |>
    subset(select = c(sampID, ct, ct_sub, taxa_id, biomass, biomass_sub,
                      density, is_distinct, is_distinct_sub,
                      'SUBORDER', 'SUBCLASS', 'PHYLUM'))

  divMets.full.rot <- calcZoopDivMetrics(subset(indata.taxa, PHYLUM=='ROTIFERA'),
                                         sampID, is_distinct, ct, suffix = 'ROT')

  divMets.full.clad <- calcZoopDivMetrics(subset(indata.taxa, SUBORDER=='CLADOCERA'),
                                         sampID, is_distinct, ct, suffix = 'CLAD')

  divMets.full.cope <- calcZoopDivMetrics(subset(indata.taxa, SUBCLASS=='COPEPODA'),
                                         sampID, is_distinct, ct, suffix = 'COPE')

  divMets.sub.rot <- calcZoopDivMetrics(subset(indata.taxa, PHYLUM=='ROTIFERA'),
                                         sampID, is_distinct_sub, ct_sub, suffix = 'ROT')

  divMets.sub.rot$PARAMETER <- paste0(divMets.sub.rot$PARAMETER, sub_mod)

  divMets.sub.clad <- calcZoopDivMetrics(subset(indata.taxa, SUBORDER=='CLADOCERA'),
                                          sampID, is_distinct_sub, ct_sub, suffix = 'CLAD')

  divMets.sub.clad$PARAMETER <- paste0(divMets.sub.clad$PARAMETER, sub_mod)

  divMets.sub.cope <- calcZoopDivMetrics(subset(indata.taxa, SUBCLASS=='COPEPODA'),
                                          sampID, is_distinct_sub, ct_sub, suffix = 'COPE')

  divMets.sub.cope$PARAMETER <- paste0(divMets.sub.cope$PARAMETER, sub_mod)

  print("Finished group diversity metrics")
# Dominance metrics
  dom.full <- calcZoopDomMetrics(indata, sampID, is_distinct,
                                  valsIn = c(ct, biomass, density),
                                  valsOut = c('PIND', 'PBIO', 'PDEN'),
                                  taxa_id,
                                  subgrp = NULL)

  dom.sub <- calcZoopDomMetrics(indata, sampID, is_distinct_sub,
                                 valsIn = c(ct_sub, biomass_sub),
                                 valsOut = c(paste(sub_mod, 'PIND', sep = '_'),
                                             paste(sub_mod, 'PBIO', sep = '_'),
                                             paste(sub_mod, 'PDEN', sep = '_')),
                                 taxa_id,
                                 subgrp = NULL)

  indata.taxa$ROT <- with(indata.taxa, ifelse(PHYLUM=='ROTIFERA', 1, 0))
  indata.taxa$CLAD <- with(indata.taxa, ifelse(SUBORDER=='CLADOCERA', 1, 0))
  indata.taxa$COPE <- with(indata.taxa, ifelse(SUBCLASS=='COPEPODA', 1, 0))

  dom.full.rot <- calcZoopDomMetrics(indata.taxa,
                                     sampID, is_distinct,
                                     valsIn = c(ct, biomass, density),
                                     valsOut = c('PIND', 'PBIO', 'PDEN'),
                                     taxa_id, subgrp = 'ROT')

  dom.full.clad <- calcZoopDomMetrics(indata.taxa,
                                     sampID, is_distinct,
                                     valsIn = c(ct, biomass, density),
                                     valsOut = c('PIND', 'PBIO', 'PDEN'),
                                     taxa_id, subgrp = 'CLAD')

  dom.full.cope <- calcZoopDomMetrics(indata.taxa,
                                     sampID, is_distinct,
                                     valsIn = c(ct, biomass, density),
                                     valsOut = c('PIND', 'PBIO', 'PDEN'),
                                     taxa_id, subgrp = 'COPE')

  dom.sub.rot <- calcZoopDomMetrics(indata.taxa,
                                     sampID, is_distinct_sub,
                                     valsIn = c(ct_sub, biomass_sub),
                                     valsOut = c('PIND', 'PBIO'),
                                     taxa_id, subgrp = 'ROT')
  dom.sub.rot$PARAMETER <- gsub('ROT\\_PIND', paste(sub_mod, 'ROT\\_PIND', sep = '_'),
                                dom.sub.rot$PARAMETER)
  dom.sub.rot$PARAMETER <- gsub('ROT\\_PBIO', paste(sub_mod, 'ROT\\_PBIO', sep = '_'),
                                dom.sub.rot$PARAMETER)

  dom.sub.clad <- calcZoopDomMetrics(subset(indata.taxa, SUBORDER=='CLADOCERA'),
                                     sampID, is_distinct_sub,
                                     valsIn = c(ct_sub, biomass_sub),
                                      valsOut = c('PIND', 'PBIO'),
                                      taxa_id, subgrp = 'CLAD')
  dom.sub.clad$PARAMETER <- gsub('CLAD\\_PIND', paste(sub_mod, 'CLAD\\_PIND', sep = '_'),
                                 dom.sub.clad$PARAMETER)
  dom.sub.clad$PARAMETER <- gsub('CLAD\\_PBIO', paste(sub_mod, 'CLAD\\_PBIO', sep = '_'),
                                 dom.sub.clad$PARAMETER)

  dom.sub.cope <- calcZoopDomMetrics(subset(indata.taxa, SUBCLASS=='COPEPODA'),
                                     sampID, is_distinct_sub,
                                     valsIn = c(ct_sub, biomass_sub),
                                     valsOut = c('PIND', 'PBIO'),
                                     taxa_id, subgrp = 'COPE')
  dom.sub.cope$PARAMETER <- gsub('COPE\\_PIND', paste(sub_mod, 'COPE\\_PIND', sep = '_'),
                                 dom.sub.cope$PARAMETER)
  dom.sub.cope$PARAMETER <- gsub('COPE\\_PBIO', paste(sub_mod, 'COPE\\_PBIO', sep = '_'),
                                 dom.sub.cope$PARAMETER)

  print("Finished dominance metrics")

  # Richness metrics
  richMets.full <- calcZoopRichnessMetrics(indata, sampID,
                                           distVars=c(is_distinct, is_distinct_sub),
                                           nonnative, inTaxa,
                                           taxa_id, genus,
                                           family, prefix = c('', sub_mod))

  print("Finished richness metrics")

  # Totals - these are all in wide format so combine and then melt
  totMets.full <- calcZoopTotals(indata, sampID, is_distinct,
                                 c(ct, biomass, density),
                                 c('TOTL_NIND', 'TOTL_BIO', 'TOTL_DEN'),
                                 'TOTL_NTAX')

  totMets.nat <- calcZoopTotals(indata.nat, sampID, is_distinct,
                                c(ct, biomass, density),
                                c('TOTL_NAT_NIND', 'TOTL_NAT_BIO', 'TOTL_NAT_DEN'),
                                'TOTL_NAT_NTAX')

  totMets.sub <- calcZoopTotals(indata, sampID, is_distinct_sub,
                                 c(ct_sub, biomass_sub),
                                 c(paste0('TOTL', sub_mod, '_NIND'),
                                   paste0('TOTL', sub_mod, '_BIO')),
                                 paste0('TOTL', sub_mod, '_NTAX'))

  totMets.sub.nat <- calcZoopTotals(indata.nat, sampID, is_distinct_sub,
                                c(ct_sub, biomass_sub),
                                c(paste0('TOTL', sub_mod, '_NAT_NIND'),
                                  paste0('TOTL', sub_mod, '_NAT_BIO')),
                                  paste0('TOTL', sub_mod, '_NAT_NTAX'))

  print("Finished totals metrics for combined sample")

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

  totMets.zocn.sub <- calcZoopTotals(inCoarse, sampID, is_distinct_sub,
                                     c(ct_sub, biomass_sub),
                                     c(paste0('ZOCN', sub_mod, '_NIND'),
                                       paste0('ZOCN', sub_mod, '_BIO')),
                                       paste0('ZOCN', sub_mod, '_NTAX'))

  totMets.zocn.sub.nat <- calcZoopTotals(subset(inCoarse, eval(as.name(nonnative))==0),
                                         sampID, is_distinct_sub,
                                         c(ct_sub, biomass_sub),
                                         c(paste0('ZOCN', sub_mod, '_NAT_NIND'),
                                           paste0('ZOCN', sub_mod, '_NAT_BIO')),
                                           paste0('ZOCN', sub_mod, '_NAT_NTAX'))

  totMets.zofn.full <- calcZoopTotals(inFine, sampID, is_distinct,
                                      c(ct, biomass, density),
                                      c('ZOFN_NIND', 'ZOFN_BIO', 'ZOFN_DEN'),
                                      'ZOFN_NTAX')

  totMets.zofn.nat <- calcZoopTotals(subset(inFine, eval(as.name(nonnative))==0),
                                     sampID, is_distinct,
                                     c(ct, biomass, density),
                                     c('ZOFN_NAT_NIND', 'ZOFN_NAT_BIO', 'ZOFN_NAT_DEN'),
                                     'ZOFN_NAT_NTAX')

  totMets.zofn.sub <- calcZoopTotals(inFine, sampID, is_distinct_sub,
                                     c(ct_sub, biomass_sub),
                                     c(paste0('ZOFN', sub_mod, '_NIND'),
                                       paste0('ZOFN', sub_mod, '_BIO')),
                                       paste0('ZOFN', sub_mod, '_NTAX'))

  totMets.zofn.sub.nat <- calcZoopTotals(subset(inFine, eval(as.name(nonnative))==0),
                                         sampID, is_distinct_sub,
                                         c(ct_sub, biomass_sub),
                                         c(paste0('ZOFN', sub_mod, '_NAT_NIND'),
                                           paste0('ZOFN', sub_mod, '_NAT_BIO')),
                                           paste0('ZOFN', sub_mod, '_NAT_NTAX'))

  print("Finished totals metrics for fine and coarse mesh data")

  zonwMets.all <- merge(totMets.full, totMets.nat, by = sampID, all.x=TRUE) |>
    merge(totMets.sub, by = sampID, all.x=TRUE) |>
    merge(totMets.sub.nat, by = sampID, all.x=TRUE)

  zocnTot.all <- merge(totMets.zocn.full, totMets.zocn.nat, by = sampID, all.x=TRUE) |>
    merge(totMets.zocn.sub, by = sampID, all.x=TRUE) |>
    merge(totMets.zocn.sub.nat, by = sampID, all.x=TRUE)

  zofnTot.all <- merge(totMets.zofn.full, totMets.zofn.nat, by = sampID, all.x=TRUE) |>
    merge(totMets.zofn.sub, by = sampID, all.x=TRUE) |>
    merge(totMets.zofn.sub.nat, by = sampID, all.x=TRUE)

  totMets.all <- merge(zonwMets.all, zocnTot.all, by = sampID, all.x=TRUE) |>
    merge(zofnTot.all, by = sampID, all.x=TRUE)

  totMets.all.long <- reshape(totMets.all, idvar = sampID, direction = 'long',
            varying = names(totMets.all)[!(names(totMets.all) %in% sampID)],
            timevar = 'PARAMETER', v.names = 'RESULT',
            times = names(totMets.all)[!(names(totMets.all) %in% sampID)])

  print("Combined totals metrics")

  # Native metrics
  natMets.full <- calcZoopNativeMetrics(totMets.all, c('UID', 'SAMPLE_TYPE'),
                        inputNative = c('TOTL_NAT_NTAX',
                                        paste0('TOTL', sub_mod, '_NAT_NTAX'),
                                        'ZOCN_NAT_NTAX',
                                        paste0('ZOCN', sub_mod, '_NAT_NTAX'),
                                        'ZOFN_NAT_NTAX',
                                        paste0('ZOFN', sub_mod, '_NAT_NTAX'),
                                        'TOTL_NAT_DEN', 'TOTL_NAT_BIO',
                                        'TOTL_NAT_NIND',
                                        paste0('TOTL', sub_mod, '_NAT_NIND'),
                                        paste0('TOTL', sub_mod, '_NAT_BIO'),
                                        'ZOCN_NAT_NIND',
                                        paste0('ZOCN', sub_mod, '_NAT_NIND'),
                                        'ZOCN_NAT_BIO',
                                        paste0('ZOCN', sub_mod, '_NAT_BIO'),
                                        'ZOCN_NAT_DEN',
                                        'ZOFN_NAT_NIND',
                                        paste0('ZOFN', sub_mod, '_NAT_NIND'),
                                        'ZOFN_NAT_BIO',
                                        paste0('ZOFN', sub_mod, '_NAT_BIO')),
                        inputTotals = c('TOTL_NTAX', paste0('TOTL', sub_mod, '_NTAX'),
                                        'ZOCN_NTAX', paste0('ZOCN', sub_mod, '_NTAX'),
                                        'ZOFN_NTAX', paste0('ZOFN', sub_mod, '_NTAX'),
                                        'TOTL_DEN',
                                        'TOTL_BIO', 'TOTL_NIND',
                                        paste0('TOTL', sub_mod, '_NIND'),
                                        paste0('TOTL', sub_mod, '_BIO'),
                                        'ZOCN_NIND',
                                        paste0('ZOCN', sub_mod, '_NIND'),
                                        'ZOCN_BIO',
                                        paste0('ZOCN', sub_mod, '_BIO'),
                                        'ZOCN_DEN', 'ZOFN_NIND',
                                        paste0('ZOFN', sub_mod, '_NIND'),
                                        'ZOFN_BIO',
                                        paste0('ZOFN', sub_mod, '_BIO')))

  print("Finished native metrics")

  # Now combine all metrics together
  allMets <- rbind(baseMets.full, baseMets.nat, baseMets.sub,
                   baseMets.sub.nat, copeRatMets,
                   divMets.full, divMets.sub,
                   divMets.full.rot, divMets.full.clad,
                   divMets.full.cope, divMets.sub.rot,
                   divMets.sub.clad, divMets.sub.cope,
                   dom.full, dom.sub, dom.full.rot,
                   dom.full.clad, dom.full.cope,
                   dom.sub.rot, dom.sub.clad,
                   dom.sub.cope, richMets.full,
                   totMets.all.long, natMets.full)
}
