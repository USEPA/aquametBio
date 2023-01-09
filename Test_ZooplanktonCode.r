# Test_ZooplanktonCode.r
# Use this code to pull data and test functions to make sure the outputs
# match up with values in database.
options(scipen = 999)

library(reshape2)
library(plyr)
library(dplyr)
library(gtools)
library(Hmisc)
library(lubridate)
library(stringr)
library(aquametBio)
library(tidyr)
library(RODBC)

source( 'l:/Priv/CORFiles/IM/Rwork/SharedCode/db.r')
source( 'l:/Priv/CORFiles/IM/Rwork/SharedCode/dd.r')
source( 'l:/Priv/CORFiles/IM/Rwork/SharedCode/md.r')
source( 'l:/Priv/CORFiles/IM/Rwork/SharedCode/sharedSupport.r')
source("R/rarifyCounts.r")
source("R/convertZoopCts_NLA.r")
source("R/prepZoopCombCts_NLA.r")
source("R/calcZoopTotals.r")
source("R/calcZoopBaseMetrics.r")
source("R/calcZoopDivMetrics.r")

rawZp <- dbGet('ALL_THE_NLA', 'tblZOOPRAW', where = "PARAMETER IN('ABUNDANCE_TOTAL', 'BIOMASS_FACTOR', 'CONCENTRATED_VOLUME', 'VOLUME_COUNTED', 'L_R_ABUND', 'LARGE_RARE_TAXA') AND
               SAMPLE_TYPE IN('ZOCN', 'ZOFN')") %>%
  pivot_wider(id_cols = c('UID', 'SAMPLE_TYPE', 'TAXA_ID', 'REP'), names_from='PARAMETER', values_from='RESULT')

towvol <- dbGet('ALL_THE_NLA', 'tblSAMPLES', where = "SAMPLE_TYPE IN('ZOCN', 'ZOFN') AND PARAMETER='TOW_VOLUME'") %>%
  pivot_wider(id_cols = c('UID', 'SAMPLE_TYPE'), names_from='PARAMETER', values_from='RESULT')

rawZp.1 <- merge(rawZp, towvol, by=c('UID', 'SAMPLE_TYPE'))
# Convert raw counts (which may have more than one rep if same taxon in count more than once)
zpCts <- convertZoopCts_NLA(rawZp.1, sampID='UID', sampType='SAMPLE_TYPE',
                            rawCt = 'ABUNDANCE_TOTAL', biofactor = 'BIOMASS_FACTOR',
                            tow_vol = 'TOW_VOLUME', vol_ctd = 'VOLUME_COUNTED',
                            conc_vol = 'CONCENTRATED_VOLUME',
                            taxa_id = 'TAXA_ID', subsample=FALSE)
# This version includes Large-rare taxa in the data (should only be for ZOCN samples)
zpCts.lr <- convertZoopCts_NLA(rawZp.1, sampID='UID', sampType='SAMPLE_TYPE',
                            rawCt = 'ABUNDANCE_TOTAL', biofactor = 'BIOMASS_FACTOR',
                            tow_vol = 'TOW_VOLUME', vol_ctd = 'VOLUME_COUNTED',
                            conc_vol = 'CONCENTRATED_VOLUME', lr_taxon = 'LARGE_RARE_TAXA',
                            taxa_id = 'TAXA_ID', subsample=FALSE)


# Compare to values in database in tblZOOPCNT
curCts <- dbGet('ALL_THE_NLA', 'tblZOOPCNT', where = "PARAMETER IN('BIOMASS', 'COUNT', 'DENSITY', 'LARGE_RARE_TAXA') AND SAMPLE_TYPE IN('ZOCN', 'ZOFN')")

compareCts <- mutate(zpCts.lr, COUNT=as.character(COUNT), BIOMASS = as.character(BIOMASS),
                     DENSITY = as.character(DENSITY)) %>%
  pivot_longer(cols=COUNT:LARGE_RARE_TAXA, names_to = 'PARAMETER',
                           values_to='RESULT', values_drop_na = TRUE) %>%
  merge(curCts, by=c('UID', 'SAMPLE_TYPE', 'TAXA_ID', 'PARAMETER'), all=TRUE)

diffCts <- filter(compareCts, PARAMETER %in% c('BIOMASS', 'COUNTS', 'DENSITY')) %>%
  mutate(RESULT.x = as.numeric(RESULT.x), RESULT.y = as.numeric(RESULT.y)) %>%
  filter(abs(RESULT.x-RESULT.y)>0.01) %>%
  mutate(diffFactor = RESULT.x/RESULT.y) # 0 differences now

msgCur <- filter(compareCts, is.na(RESULT.y) & !is.na(RESULT.x))
msgNew <- filter(compareCts, is.na(RESULT.x) & !is.na(RESULT.y))

# Keep only counts where both ZOFN and ZOCN samples are in counts
numsamps <- dbGet('ALL_THE_NLA', 'tblZOOPCNT', where = "SAMPLE_TYPE IN('ZOCN', 'ZOFN')") %>%
  select(UID, SAMPLE_TYPE) %>%
  unique() %>%
  ddply(c('UID'), summarise, NUM = length(SAMPLE_TYPE))

msgSamps <- filter(numsamps, NUM == 1)

keepSamps <- filter(numsamps, NUM==2)

# Now look at code to combine values from ZOCN and ZOFN samples
zpCts.in <- filter(zpCts.lr, UID %in% keepSamps$UID)

zonwCts <- prepZoopCombCts_NLA(zpCts.in, sampID = 'UID', sampType = 'SAMPLE_TYPE', typeFine = 'ZOFN', typeCoarse = 'ZOCN', ct = 'COUNT', biomass = 'BIOMASS', density = 'DENSITY', taxa_id = 'TAXA_ID', lr_taxon = 'LARGE_RARE_TAXA')

curZONW <- dbGet('ALL_THE_NLA', 'tblZOOPCNT', where = "SAMPLE_TYPE = 'ZONW' AND PARAMETER IN('BIOMASS', 'COUNT', 'DENSITY', 'LARGE_RARE_TAXA')")

compareZONW <- mutate(zonwCts, COUNT=as.character(COUNT), BIOMASS = as.character(BIOMASS),
                      DENSITY = as.character(DENSITY)) %>%
  pivot_longer(cols = COUNT:LARGE_RARE_TAXA, names_to='PARAMETER',
                            values_to = 'RESULT', values_drop_na = TRUE) %>%
  merge(curZONW, by=c('UID', 'SAMPLE_TYPE', 'TAXA_ID', 'PARAMETER'), all=TRUE)

diffZONW <- filter(compareZONW, PARAMETER %in% c('BIOMASS', 'COUNTS', 'DENSITY')) %>%
  mutate(RESULT.x = as.numeric(RESULT.x), RESULT.y = as.numeric(RESULT.y)) %>%
  filter(abs(RESULT.x-RESULT.y)>0.01) %>%
  mutate(diffFactor = RESULT.x/RESULT.y) # 0 differences now

msgCur.zonw <- filter(compareZONW, is.na(RESULT.y)) # Only missing count for taxa_ID 9999, does not really matter because 9999 is for sample with 0 individuals in it.
msgNew.zonw <- filter(compareZONW, is.na(RESULT.x))

# Now we need to create the 300 count subsamples for each sample type, then combine only for BIOMASS and COUNT. Start with raw data
rawZoop.300 <- rarifyCounts(rawZp.1, sampID = c('UID', 'SAMPLE_TYPE'),
                            abund = 'ABUNDANCE_TOTAL',
                            subsize = 300,
                            seed = 1234)

rawZoop.300.sub <- filter(rawZoop.300, ABUNDANCE_TOTAL>0)

zpCts.300 <- convertZoopCts_NLA(rawZoop.300.sub, sampID='UID', sampType='SAMPLE_TYPE',
                                rawCt = 'ABUNDANCE_TOTAL', biofactor = 'BIOMASS_FACTOR',
                                tow_vol = 'TOW_VOLUME', vol_ctd = 'VOLUME_COUNTED',
                                conc_vol = 'CONCENTRATED_VOLUME',
                                taxa_id = 'TAXA_ID', subsample=TRUE)

zpCts.300 <- plyr::rename(zpCts.300, c('COUNT'='COUNT_300', 'BIOMASS' = 'BIOMASS_300'))

zpCts.in.300 <- filter(zpCts.300, UID %in% keepSamps$UID)

zonwCts.300 <- prepZoopCombCts_NLA(zpCts.in.300, sampID = 'UID',
                                   sampType = 'SAMPLE_TYPE', typeFine = 'ZOFN',
                                   typeCoarse = 'ZOCN', ct = 'COUNT_300',
                                   biomass = 'BIOMASS_300', taxa_id = 'TAXA_ID')

# Now determine distinctness for each separate set of counts by UID, sample type
# Combine ZONW, ZOCN, and ZOFN together for each sample size
taxa <- dbGet('ALL_THE_NLA', 'tblTAXA', where = "ASSEMBLAGE_NAME='ZOOPLANKTON'") %>%
  pivot_wider(id_cols = c('TAXA_ID'), names_from='PARAMETER', values_from='RESULT')

zpCts.all <- rbind(zonwCts, zpCts.lr) %>%
  filter(COUNT>0|!is.na(LARGE_RARE_TAXA)) %>%
  merge(taxa, by = c('TAXA_ID'))

# Add in large/rare taxa for ZOCN samples and ZONW samples and calculate distinctness, including L/R taxa

zpCts.all.dist <- assignDistinct(zpCts.all, sampleID = c('UID', 'SAMPLE_TYPE'), taxlevels=c('PHYLUM','CLASS','SUBCLASS','ORDER','SUBORDER','FAMILY','GENUS','SPECIES','SUBSPECIES'))

zpCts.all.300 <- rbind(zonwCts.300, zpCts.300) %>%
  filter(COUNT_300>0) %>%
  merge(taxa, by = c('TAXA_ID'))

zpCts.all.300.dist <- assignDistinct(zpCts.all.300, sampleID = c('UID', 'SAMPLE_TYPE'), taxlevels=c('PHYLUM','CLASS','SUBCLASS','ORDER','SUBORDER','FAMILY','GENUS','SPECIES','SUBSPECIES')) %>%
  plyr::rename(c('IS_DISTINCT' = 'IS_DISTINCT_300'))


# Compare IS_DISTINCT for full counts only (since we can't repeat the random subsampling exactly to match the database)
curDist <- dbGet('ALL_THE_NLA', 'tblZOOPCNT', where = "PARAMETER IN('COUNT', 'IS_DISTINCT')") %>%
  pivot_wider(id_cols = c('UID', 'SAMPLE_TYPE', 'TAXA_ID'), names_from='PARAMETER', values_from='RESULT')

matchDist.full <- merge(curDist, zpCts.all.dist, by=c('UID', 'SAMPLE_TYPE', 'TAXA_ID'))

diffs.Dist.full <- filter(matchDist.full, IS_DISTINCT.x!=IS_DISTINCT.y) # All of these were assigned distinctness by Dave Peck in 2012 and all have subspecies in the same sample - see below -

table(diffs.Dist.full$TAXA_ID)

curDist.sub1245 <- filter(curDist, TAXA_ID %in% c(1246, 1247, 1248, 1249)) %>%
  merge(subset(diffs.Dist.full, TAXA_ID==1245), by=c('UID', 'SAMPLE_TYPE'))

nrow(unique(curDist.sub1245[, c('UID', 'SAMPLE_TYPE')]))

curDist.sub1254 <- filter(curDist, TAXA_ID %in% c(1255)) %>%
  merge(subset(diffs.Dist.full, TAXA_ID==1254), by=c('UID', 'SAMPLE_TYPE'))

nrow(unique(curDist.sub1254[, c('UID', 'SAMPLE_TYPE')]))

curDist.sub1072 <- filter(curDist, TAXA_ID %in% c(1073, 1074)) %>%
  merge(subset(diffs.Dist.full, TAXA_ID==1072), by=c('UID', 'SAMPLE_TYPE'))

nrow(unique(curDist.sub1072[, c('UID', 'SAMPLE_TYPE')]))

verif <- dbGet('ALL_THE_NLA', 'tblVERIFICATION', where = "PARAMETER = 'DATE_COL'") %>%
  filter(UID %in% diffs.Dist.full$UID)
range(verif$RESULT)

# Cannot test for 300 subsample because we don't have exactly the same data, but we can pull the data
# from the database, re-assign IS_DISTINCT, and compare to original.
cur300 <- dbGet('ALL_THE_NLA', 'tblZOOPCNT', where = "PARAMETER IN('COUNT_300', 'IS_DISTINCT_300')")

cur300.ct <- subset(cur300, PARAMETER == 'COUNT_300') %>%
  mutate(RESULT = as.numeric(RESULT)) %>%
  filter(RESULT >0) %>%
  merge(taxa, by='TAXA_ID')

cur300.addDist <- assignDistinct(cur300.ct, sampleID = c('UID', 'SAMPLE_TYPE'), taxlevels=c('PHYLUM','CLASS','SUBCLASS','ORDER','SUBORDER','FAMILY','GENUS','SPECIES','SUBSPECIES')) %>%
  plyr::rename(c('IS_DISTINCT' = 'IS_DISTINCT_300'))

matchDist.300 <- merge(filter(cur300, PARAMETER=='IS_DISTINCT_300'), cur300.addDist, by=c('UID', 'SAMPLE_TYPE', 'TAXA_ID'))

diffs.300 <- filter(matchDist.300, RESULT.x!=IS_DISTINCT_300) # 242 records, same taxa as above (1072, 1245, 1254)


# Test totals code with various inputs and arguments
curZp <- dbGet('ALL_THE_NLA', 'tblZOOPCNT', where = "SAMPLE_TYPE IN('ZONW')") %>%
  pivot_wider(id_cols = c("UID", "SAMPLE_TYPE", "TAXA_ID"), names_from='PARAMETER', values_from='RESULT')

testSum <- calcZoopTotals(indata = curZp, sampID = c('UID', 'SAMPLE_TYPE'),
                          is_distinct = 'IS_DISTINCT',
                          inputSums = c('COUNT', 'BIOMASS', 'DENSITY'),
                          outputSums = c('TOTL_NIND', 'TOTL_BIO', 'TOTL_DEN'),
                          outputTaxa = 'TOTL_NTAX')

testSum.300 <- calcZoopTotals(curZp, c('UID', 'SAMPLE_TYPE'),
                              'IS_DISTINCT',
                            c('COUNT_300', 'BIOMASS_300'),
                            c('TOTL300_NIND', 'TOTL300_BIO'),
                            'TOTL300_NTAX')

testSum.nat <- calcZoopTotals(indata = subset(curZp, is.na(NON_NATIVE)|NON_NATIVE=='N'),
                              sampID = c('UID', 'SAMPLE_TYPE'),
                              is_distinct = 'IS_DISTINCT',
                              inputSums = c('COUNT', 'BIOMASS', 'DENSITY'),
                              outputSums = c('TOTL_NAT_NIND', 'TOTL_NAT_BIO', 'TOTL_NAT_DEN'),
                              outputTaxa = 'TOTL_NAT_NTAX')

# Run metrics code using current zp count data, then apply to native only
taxa.clean <- filter(taxa, is.na(NON_TARGET))

zpIn <- filter(curZp, TAXA_ID %in% taxa.clean$TAXA_ID)

testMets <- calcZoopMetrics(zpIn, c('UID', 'SAMPLE_TYPE'), 'IS_DISTINCT',
                                        'COUNT', 'BIOMASS', 'DENSITY',
                                        taxa, taxa_id='TAXA_ID',
                                        'FFG', 'CLADOCERA_SIZE',
                                        'NET_SIZECLS_NEW',
                                        nativeMetrics = FALSE)

# Now try with only native data
state <- dbGet('ALL_THE_NLA', 'tblSITETRUTH', where = "IND_DOMAIN IN('HAND', 'CORE') AND STUDY='NLA' AND PARAMETER IN('SITE_ID', 'PSTL_CODE')") %>%
  pivot_wider(id_cols = c('UNIQUE_ID', 'DSGN_CYCLE', 'IND_DOMAIN'), names_from='PARAMETER', values_from='RESULT')

verif <- dbGet('ALL_THE_NLA', 'tblVERIFICATION', where = "PARAMETER IN('SITE_ID')") %>%
  pivot_wider(id_cols = 'UID', names_from='PARAMETER', values_from = 'RESULT')

stateVerif <- merge(verif, state, by='SITE_ID')

zpIn.nat <- merge(curZp, taxa.clean, by='TAXA_ID') %>%
  merge(stateVerif, by=c('UID')) %>%
  filter(is.na(NON_NATIVE.x)|NON_NATIVE.x=='N') %>%
  filter(NON_NATIVE.y!= "ALL"|is.na(NON_NATIVE.y)) %>%
  filter(str_detect(NON_NATIVE.y, PSTL_CODE)==FALSE|is.na(NON_NATIVE.y)) %>%
  select(UID, TAXA_ID, SAMPLE_TYPE, COUNT, BIOMASS, DENSITY, IS_DISTINCT,
         LARGE_RARE_TAXA)

testMets.nat <- calcZoopMetrics(zpIn.nat, c('UID', 'SAMPLE_TYPE'),
                                'IS_DISTINCT',
                                'COUNT', 'BIOMASS', 'DENSITY',
                                taxa, taxa_id='TAXA_ID',
                                'FFG', 'CLADOCERA_SIZE',
                                'NET_SIZECLS_NEW',
                                nativeMetrics = TRUE)

# Compare values to database, assuming all preparatory steps are being followed
compMets <- dbGet('ALL_THE_NLA', 'tblZOOPMET') %>%
  merge(testMets, by=c('UID', 'PARAMETER'))

diffMets <- filter(compMets, abs(as.numeric(RESULT.x)-RESULT.y)>0.0001) %>%
  mutate(DIFF = as.numeric(RESULT.x) - RESULT.y) # All within rounding limits


# Zooplankton diversity metrics
divMets <- calcZoopDivMetrics(zpIn, c('UID', 'SAMPLE_TYPE'), 'IS_DISTINCT',
                              'COUNT', 'BIOMASS', 'DENSITY')

divMets.long <- pivot_longer(divMets, cols = HPRIME_NIND:SIMPSON_BIO, names_to='PARAMETER',
                             values_to='RESULT')

compDivMets <- dbGet('ALL_THE_NLA', 'tblZOOPMET') %>%
  merge(divMets.long, by=c('UID', 'PARAMETER'))

diffMets.div <- filter(compDivMets, abs(as.numeric(RESULT.x)-RESULT.y)>0.0002) %>%
  mutate(DIFF = as.numeric(RESULT.x) - RESULT.y) %>% # Only HPRIME and SIMPSON for BIO and DEN variations, up to 0.0010, these do not differ using the code
  merge(stateVerif, by='UID') # Rounding changed with version 4 of R, so that might account for some of the variation that I cannot reproduce using the same code used for the original metric values.

# indf.sums <- plyr::ddply(indata.1, c('UID'), mutate,
#                    SUMBIO=sum(as.numeric(BIOMASS), na.rm=T),
#                    SUMDEN=sum(as.numeric(DENSITY), na.rm=T))
# indf.props <- mutate(indf.sums, prop.bio = BIOMASS/SUMBIO, prop.den = DENSITY/SUMDEN)
#
# test.props <- select(indf.props, UID, TAXA_ID, prop.bio, prop.den) %>%
#   merge(select(calcData, UID, TAXA_ID, prop.bio, prop.den), by = c('UID', 'TAXA_ID'))
#
# test.div <- ddply(indf.props, c('UID'), summarise,
#                   HPRIME_BIO=round(-1*sum(prop.bio*log(prop.bio),na.rm=T),4),
#                   HPRIME_DEN=round(-1*sum(prop.den*log(prop.den),na.rm=T),4),
#                   SIMPSON_BIO=round(sum(prop.bio*prop.bio,na.rm=T),4),
#                   SIMPSON_DEN=round(sum(prop.den*prop.den,na.rm=T),4))
#
# filter(test.props, abs(prop.bio.x - prop.bio.y)>0.0001)
# filter(test.props, abs(prop.den.x - prop.den.y) >0.0001)
#
# test.divMets <- merge(divMets, test.div, by='UID')
#
# filter(test.divMets, abs(HPRIME_BIO.x-HPRIME_BIO.y)>0.0001)
# filter(test.divMets, abs(HPRIME_DEN.x-HPRIME_DEN.y)>0.0001)
# filter(test.divMets, abs(SIMPSON_BIO.x-SIMPSON_BIO.y)>0.0001)
# filter(test.divMets, abs(SIMPSON_DEN.x-SIMPSON_DEN.y)>0.0001)
