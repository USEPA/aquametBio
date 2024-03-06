## ----convertRaw-1, message=FALSE, warning=FALSE-------------------------------
library(aquametBio)

head(zpRawEx)

zpCts.full <- convertZoopCts_NLA(zpRawEx, sampID='UID', sampType='SAMPLE_TYPE',
                            rawCt = 'ABUNDANCE_TOTAL', biofactor = 'BIOMASS_FACTOR',
                            tow_vol = 'TOW_VOLUME', vol_ctd = 'VOLUME_COUNTED',
                            conc_vol = 'CONCENTRATED_VOLUME', lr_taxon = 'LARGE_RARE_TAXA',
                            taxa_id = 'TAXA_ID', subsample=FALSE)

zpCts.full <- subset(zpCts.full, COUNT!=0|LARGE_RARE_TAXA=='Y') |>
  merge(zoopTaxa, by = 'TAXA_ID') |>
  assignDistinct(sampleID = c('UID', 'SAMPLE_TYPE'), 
                 taxlevels=c('PHYLUM','CLASS','SUBCLASS',
                             'ORDER','SUBORDER','FAMILY','GENUS',
                             'SPECIES','SUBSPECIES')) 

head(zpCts.full)

## ----createZONW, message=FALSE, warning=FALSE---------------------------------
zpCts.zonw <- prepZoopCombCts_NLA(zpCts.full, sampID = c('UID'),
                                  sampType = 'SAMPLE_TYPE', typeFine = 'ZOFN',
                                  typeCoarse = 'ZOCN', ct = 'COUNT',
                                  biomass = 'BIOMASS', density = 'DENSITY', 
                                  taxa_id = 'TAXA_ID', lr_taxon = 'LARGE_RARE_TAXA')

zpCts.zonw <- subset(zpCts.zonw, COUNT!=0|LARGE_RARE_TAXA=='Y') |>
  merge(zoopTaxa, by = 'TAXA_ID') |>
  assignDistinct(sampleID = c('UID', 'SAMPLE_TYPE'), 
                 taxlevels=c('PHYLUM','CLASS','SUBCLASS',
                             'ORDER','SUBORDER','FAMILY','GENUS',
                             'SPECIES','SUBSPECIES')) 


head(zpCts.zonw)

## ----rarify, message=FALSE, warning = FALSE-----------------------------------
zpRawEx.300 <- rarifyCounts(zpRawEx, sampID = c('UID', 'SAMPLE_TYPE'),
                            abund = 'ABUNDANCE_TOTAL',
                            subsize = 300,
                            seed = 1234) |>
  subset(ABUNDANCE_TOTAL > 0, select = -LARGE_RARE_TAXA)

head(zpRawEx.300)

## ----subsamp-1, message=FALSE, warning = FALSE--------------------------------
zpCts.300 <- convertZoopCts_NLA(zpRawEx.300, sampID='UID', sampType='SAMPLE_TYPE',
                            rawCt = 'ABUNDANCE_TOTAL', biofactor = 'BIOMASS_FACTOR',
                            tow_vol = 'TOW_VOLUME', vol_ctd = 'VOLUME_COUNTED',
                            conc_vol = 'CONCENTRATED_VOLUME', 
                            taxa_id = 'TAXA_ID', subsample=TRUE) 

zpCts.300 <- subset(zpCts.300, COUNT!=0) |>
  merge(zoopTaxa, by = 'TAXA_ID') |>
  assignDistinct(sampleID = c('UID', 'SAMPLE_TYPE'), 
                 taxlevels=c('PHYLUM','CLASS','SUBCLASS',
                             'ORDER','SUBORDER','FAMILY','GENUS',
                             'SPECIES','SUBSPECIES')) |>
  subset(select = c('UID', 'SAMPLE_TYPE', 'TAXA_ID', 'COUNT', 'BIOMASS', 'IS_DISTINCT'))

names(zpCts.300)[names(zpCts.300) %in% c('COUNT', 'BIOMASS', 'IS_DISTINCT')] <- paste(names(zpCts.300)[names(zpCts.300) %in% c('COUNT', 'BIOMASS', 'IS_DISTINCT')], "300", sep="_") 

head(zpCts.300)

## ----subsamp-comb, message=FALSE, warning = FALSE-----------------------------
zpCts.zonw.300 <- prepZoopCombCts_NLA(zpCts.300, sampID = c('UID'),
                                  sampType = 'SAMPLE_TYPE', typeFine = 'ZOFN',
                                  typeCoarse = 'ZOCN', ct = 'COUNT_300',
                                  biomass = 'BIOMASS_300', 
                                  taxa_id = 'TAXA_ID')

zpCts.zonw.300 <- subset(zpCts.zonw.300, COUNT_300!=0) |>
  merge(zoopTaxa, by = 'TAXA_ID') |>
  assignDistinct(sampleID = c('UID', 'SAMPLE_TYPE'), 
                 taxlevels=c('PHYLUM','CLASS','SUBCLASS',
                             'ORDER','SUBORDER','FAMILY','GENUS',
                             'SPECIES','SUBSPECIES')) |>
  subset(select = c('UID', 'SAMPLE_TYPE', 'TAXA_ID', 'COUNT_300', 'BIOMASS_300', 'IS_DISTINCT'))

names(zpCts.zonw.300)[names(zpCts.zonw.300) %in% c('IS_DISTINCT')] <- 'IS_DISTINCT_300'

head(zpCts.zonw.300)


## ----combine, message=FALSE, warning = FALSE----------------------------------
zpCts.zonw.all <- merge(zpCts.zonw, zpCts.zonw.300, 
                   by = c('UID', 'SAMPLE_TYPE', 'TAXA_ID'), all.x=TRUE)


head(zpCts.zonw.all)

zpCts.full.zocn <- subset(zpCts.full, SAMPLE_TYPE=='ZOCN')
zpCts.300.zocn <- subset(zpCts.300, SAMPLE_TYPE=='ZOCN')

zpCts.zocn.all <- merge(zpCts.full.zocn, zpCts.300.zocn, 
                  by = c('UID', 'SAMPLE_TYPE', 'TAXA_ID'), all.x=TRUE)

head(zpCts.zocn.all)

zpCts.full.zofn <- subset(zpCts.full, SAMPLE_TYPE=='ZOFN')
zpCts.300.zofn <- subset(zpCts.300, SAMPLE_TYPE=='ZOFN')

zpCts.zofn.all <- merge(zpCts.full.zofn, zpCts.300.zofn, 
                  by = c('UID', 'SAMPLE_TYPE', 'TAXA_ID'), all.x=TRUE)

head(zpCts.zofn.all)

## ----non-native, message=FALSE, warning = FALSE-------------------------------
states <- data.frame(UID = c(7100, 7122, 7477, 7548, 7748, 8395, 8499, 
                             2010287, 2010569, 2010833),
                     STATE = c('OR', 'PA', 'ID', 'IN', 'MN', 'UT', 'IL', 'MI', 
                               'SD', 'SC'),
                     ECO_BIO = c('WMTNS', 'EHIGH', 'WMTNS', 'EHIGH', 'UMW',
                                 'WMTNS', 'PLAINS', 'UMW', 'PLAINS', 'EHIGH'))

# Combined sample
zonwIn.nn <- merge(zpCts.zonw.all, states, by = 'UID') 

zonwIn.nn$NON_NATIVE_TAXON <- with(zonwIn.nn, 
                                          ifelse(is.na(NON_NATIVE), 0,
                                            ifelse(NON_NATIVE=='ALL', 1,
                                              ifelse(grep(STATE, NON_NATIVE), 1, 0)))) 

zonwIn.nn <- zonwIn.nn[, c('UID', 'SAMPLE_TYPE', 'TAXA_ID', 
                                           'COUNT', 'BIOMASS', 'DENSITY',
                                           'IS_DISTINCT', 'LARGE_RARE_TAXA',
                                           'COUNT_300', 'BIOMASS_300',
                                           'IS_DISTINCT_300',
                                           'ECO_BIO', 'NON_NATIVE_TAXON')]

head(zonwIn.nn)

# Coarse mesh net sample
zocnIn.nn <- merge(zpCts.zocn.all, states, by = 'UID') 

zocnIn.nn$NON_NATIVE_TAXON <- with(zocnIn.nn, 
                                          ifelse(is.na(NON_NATIVE), 0,
                                            ifelse(NON_NATIVE=='ALL', 1,
                                              ifelse(grep(STATE, NON_NATIVE), 1, 0)))) 

zocnIn.nn <- zocnIn.nn[, c('UID', 'SAMPLE_TYPE', 'TAXA_ID', 
                                           'COUNT', 'BIOMASS', 'DENSITY',
                                           'IS_DISTINCT', 'LARGE_RARE_TAXA',
                                           'COUNT_300', 'BIOMASS_300',
                                           'IS_DISTINCT_300',
                                           'ECO_BIO', 'NON_NATIVE_TAXON')]

head(zocnIn.nn) 

# Fine mesh net sample
zofnIn.nn <- merge(zpCts.zofn.all, states, by = 'UID') 

zofnIn.nn$NON_NATIVE_TAXON <- with(zofnIn.nn, 
                                          ifelse(is.na(NON_NATIVE), 0,
                                            ifelse(NON_NATIVE=='ALL', 1,
                                              ifelse(grep(STATE, NON_NATIVE), 1, 0)))) 

zofnIn.nn <- zofnIn.nn[, c('UID', 'SAMPLE_TYPE', 'TAXA_ID', 
                                           'COUNT', 'BIOMASS', 'DENSITY',
                                           'IS_DISTINCT', 'LARGE_RARE_TAXA',
                                           'COUNT_300', 'BIOMASS_300',
                                           'IS_DISTINCT_300',
                                           'ECO_BIO', 'NON_NATIVE_TAXON')]

head(zofnIn.nn)

## ----calc-mets-1.1, message=FALSE, warning = FALSE----------------------------
# Full counts with native taxa only
zoopMets.base.nat <- calcZoopBaseMetrics(indata = zonwIn.nn, 
                                         sampID = c('UID', 'SAMPLE_TYPE'), 
                                         is_distinct = 'IS_DISTINCT',
                                         ct = 'COUNT', 
                                         biomass = 'BIOMASS', 
                                         density = 'DENSITY',
                                         inTaxa = zoopTaxa, 
                                         taxa_id='TAXA_ID',
                                         ffg = 'FFG', 
                                         clad_size = 'CLADOCERA_SIZE',
                                         net_size = 'NET_SIZECLS_NEW',
                                         nativeMetrics = TRUE, 
                                         nonnative = 'NON_NATIVE_TAXON')

noquote(unique(zoopMets.base.nat$PARAMETER))

## ----calc-mets-1.2, message=FALSE, warning = FALSE----------------------------
# 300-organism subsamples
zoopMets.300 <- calcZoopBaseMetrics(indata = zonwIn.nn, 
                                    sampID = c('UID', 'SAMPLE_TYPE'), 
                                    is_distinct = 'IS_DISTINCT_300',
                                    ct = 'COUNT_300', 
                                    biomass = 'BIOMASS_300',
                                    inTaxa = zoopTaxa, 
                                    taxa_id='TAXA_ID',
                                    ffg = 'FFG', 
                                    clad_size = 'CLADOCERA_SIZE',
                                    net_size = 'NET_SIZECLS_NEW',
                                    nativeMetrics = FALSE)

noquote(unique(zoopMets.base.nat$PARAMETER))

## ----calc-mets-2.1, message=FALSE, warning = FALSE----------------------------
# Full counts with native taxa only
zoopMets.div <- calcZoopDivMetrics(zonwIn.nn, 
                                        sampID = 'UID', 
                                        is_distinct = 'IS_DISTINCT',
                                        ct = 'COUNT', 
                                        biomass = 'BIOMASS', 
                                        density = 'DENSITY')

noquote(unique(zoopMets.div$PARAMETER))


## ----calc-mets-2.2, message=FALSE, warning = FALSE----------------------------
# 300-organism subsamples
zoopMets.div.300 <- calcZoopDivMetrics(zonwIn.nn, 
                                       sampID = 'UID', 
                                       is_distinct = 'IS_DISTINCT_300',
                                       ct = 'COUNT_300', 
                                       biomass = 'BIOMASS_300')

# Update parameter names to reflect subsample
zoopMets.div.300$PARAMETER <- gsub("\\_NIND", "300\\_NIND", zoopMets.div.300$PARAMETER)
zoopMets.div.300$PARAMETER <- gsub("\\_BIO", "300\\_BIO", zoopMets.div.300$PARAMETER)

noquote(unique(zoopMets.div.300$PARAMETER))

## ----calc-mets-2.3, message=FALSE, warning = FALSE----------------------------
# Rotifer diversity
# First need to merge taxa list with data to be able to select rotifers only
zonwIn.taxa <- merge(zonwIn.nn, zoopTaxa, by = 'TAXA_ID') |>
    subset(select = c('UID', 'TAXA_ID', 'COUNT', 'COUNT_300', 'BIOMASS', 
                      'BIOMASS_300', 'DENSITY', 'IS_DISTINCT', 
                      'IS_DISTINCT_300', 'SUBORDER', 'SUBCLASS', 'PHYLUM'))

# Calculate rotifer diversity for full counts
zoopMets.div.rot <- calcZoopDivMetrics(subset(zonwIn.taxa, PHYLUM=='ROTIFERA'),
                                       sampID = 'UID', 
                                       is_distinct = 'IS_DISTINCT',
                                       ct = 'COUNT', 
                                       suffix = 'ROT')

noquote(unique(zoopMets.div.rot$PARAMETER))

## ----calc-mets-2.4, message=FALSE, warning = FALSE----------------------------
# Calculate rotifer diversity for subsample data
zoopMets.div.rot.300 <- calcZoopDivMetrics(subset(zonwIn.taxa, PHYLUM=='ROTIFERA'),
                                       sampID = 'UID', 
                                       is_distinct = 'IS_DISTINCT_300',
                                       ct = 'COUNT_300', 
                                       suffix = 'ROT')
# Update parameter names to reflect subsample
zoopMets.div.rot.300$PARAMETER <- paste0(zoopMets.div.rot.300$PARAMETER, '300')

noquote(unique(zoopMets.div.rot.300$PARAMETER))

## ----calc-mets-3.1, message = FALSE, warning = FALSE--------------------------
zoopMets.zonw <- calcZoopRichnessMetrics(zonwIn.nn, 
                                         sampID = 'UID',
                                         distVars=c('IS_DISTINCT', 'IS_DISTINCT_300'),
                                         nonnative = 'NON_NATIVE_TAXON', 
                                         inTaxa = zoopTaxa,
                                         taxa_id = 'TAXA_ID', 
                                         genus = 'GENUS',
                                         family = 'FAMILY', 
                                         prefix = c('', '300')) 

noquote(unique(zoopMets.zonw$PARAMETER))
    

## ----calc-mets-3.2, message = FALSE, warning = FALSE--------------------------
  zoopMets.zocn <- calcZoopRichnessMetrics(zocnIn.nn, 
                                             sampID = 'UID',
                                             distVars=c('IS_DISTINCT', 'IS_DISTINCT_300'),
                                             nonnative = 'NON_NATIVE_TAXON', 
                                             inTaxa = zoopTaxa,
                                             taxa_id = 'TAXA_ID', 
                                             genus = 'GENUS',
                                             family = 'FAMILY', 
                                             prefix = c('', '300'))

  zoopMets.zocn$PARAMETER <- gsub('FAM\\_', 'ZOCN\\_FAM\\_', zoopMets.zocn$PARAMETER)
  zoopMets.zocn$PARAMETER <- gsub('GEN\\_', 'ZOCN\\_GEN\\_', zoopMets.zocn$PARAMETER)
  zoopMets.zocn$PARAMETER <- gsub('FAM300', 'ZOCN300\\_FAM',
                                    zoopMets.zocn$PARAMETER)
  zoopMets.zocn$PARAMETER <- gsub('GEN300', 'ZOCN300\\_GEN',
                                    zoopMets.zocn$PARAMETER)
  
  noquote(unique(zoopMets.zocn$PARAMETER))

## ----calc-mets-all, message=FALSE, warning=FALSE------------------------------
zoopMets <- calcZoopAllMets(indata = zonwIn.nn,
                               inCoarse = zocnIn.nn,
                               inFine = zofnIn.nn,
                               inTaxa = zoopTaxa,
                               sampID = c('UID'),
                               is_distinct = 'IS_DISTINCT',
                               ct = 'COUNT', biomass = 'BIOMASS',
                               density = 'DENSITY',
                               is_distinct_sub = 'IS_DISTINCT_300',
                               ct_sub = 'COUNT_300',
                               biomass_sub = 'BIOMASS_300',
                               sub_mod = '300', taxa_id = 'TAXA_ID',
                               nonnative = 'NON_NATIVE_TAXON', genus = 'GENUS',
                               family = 'FAMILY', ffg = 'FFG',
                               clad_size = 'CLADOCERA_SIZE',
                               net_size = 'NET_SIZECLS_NEW')

head(zoopMets)

noquote(unique(zoopMets$PARAMETER))


## ----calc-MMI, message = FALSE, warning = FALSE-------------------------------
mmiIn <- merge(states, zoopMets, by = 'UID') |>
  reshape(idvar = c('UID', 'ECO_BIO'), direction = 'wide',
                             timevar = 'PARAMETER', v.names = 'RESULT')

names(mmiIn) <- gsub('RESULT\\.', '', names(mmiIn))

mmiOut <- calcNLA_ZoopMMI(inMets = mmiIn, 
                          sampID = 'UID',
                          ecoreg = 'ECO_BIO',
                          totlnind = 'TOTL_NIND')

head(mmiOut)


