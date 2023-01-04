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
curZp <- dbGet('ALL_THE_NLA', 'tblZOOPCNT', where = "SAMPLE_TYPE IN('ZOCN')") %>%
  pivot_wider(id_cols = c("UID", "SAMPLE_TYPE", "TAXA_ID"), names_from='PARAMETER', values_from='RESULT')

source("R/calcZoopTotals.r")

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
