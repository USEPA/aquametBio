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
  pivot_wider(id_cols = c('UID', 'SAMPLE_TYPE', 'TAXA_ID', 'REP'), names_from='PARAMETER', values_from='RESULT') %>%
  filter(is.na(L_R_ABUND))

towvol <- dbGet('ALL_THE_NLA', 'tblSAMPLES', where = "SAMPLE_TYPE IN('ZOCN', 'ZOFN') AND PARAMETER='TOW_VOLUME'") %>%
  pivot_wider(id_cols = c('UID', 'SAMPLE_TYPE'), names_from='PARAMETER', values_from='RESULT')

rawZp.1 <- merge(rawZp, towvol, by=c('UID', 'SAMPLE_TYPE'))
# Convert raw counts (which may have more than one rep if same taxon in count more than once)
zpCts <- convertZoopCts_NLA(rawZp.1, sampID='UID', sampType='SAMPLE_TYPE',
                            rawCt = 'ABUNDANCE_TOTAL', biofactor = 'BIOMASS_FACTOR',
                            tow_vol = 'TOW_VOLUME', vol_ctd = 'VOLUME_COUNTED',
                            conc_vol = 'CONCENTRATED_VOLUME',
                            taxa_id = 'TAXA_ID', subsample=FALSE)

# Compare to values in database in tblZOOPCNT
curCts <- dbGet('ALL_THE_NLA', 'tblZOOPCNT', where = "PARAMETER IN('BIOMASS', 'COUNT', 'DENSITY') AND SAMPLE_TYPE IN('ZOCN', 'ZOFN')")

compareCts <- pivot_longer(zpCts, cols=COUNT:DENSITY, names_to = 'PARAMETER', values_to='RESULT') %>%
  merge(curCts, by=c('UID', 'SAMPLE_TYPE', 'TAXA_ID', 'PARAMETER'), all=TRUE) %>%
  mutate(RESULT.x=as.numeric(RESULT.x), RESULT.y=as.numeric(RESULT.y))

diffCts <- filter(compareCts, abs(RESULT.x-RESULT.y)>0.01) %>%
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
zpCts.in <- filter(zpCts, UID %in% keepSamps$UID)

zonwCts <- prepZoopCombCts_NLA(zpCts.in, sampID = 'UID', sampType = 'SAMPLE_TYPE', typeFine = 'ZOFN', typeCoarse = 'ZOCN', ct = 'COUNT', biomass = 'BIOMASS', density = 'DENSITY', taxa_id = 'TAXA_ID')

curZONW <- dbGet('ALL_THE_NLA', 'tblZOOPCNT', where = "SAMPLE_TYPE = 'ZONW' AND PARAMETER IN('BIOMASS', 'COUNT', 'DENSITY')")

compareZONW <- pivot_longer(zonwCts, cols = COUNT:DENSITY, names_to='PARAMETER',
                            values_to = 'RESULT') %>%
  merge(curZONW, by=c('UID', 'SAMPLE_TYPE', 'TAXA_ID', 'PARAMETER'), all=TRUE) %>%
  mutate(RESULT.x=as.numeric(RESULT.x), RESULT.y=as.numeric(RESULT.y))

diffZONW <- filter(compareZONW, abs(RESULT.x-RESULT.y)>0.01) %>%
  mutate(diffFactor = RESULT.x/RESULT.y) # 0 differences now

msgCur.zonw <- filter(compareZONW, is.na(RESULT.y))
msgNew.zonw <- filter(compareZONW, is.na(RESULT.x)) # Only missing count for taxa_ID 9999, does not really matter because 9999 is for sample with 0 individuals in it.

# Now we need to create the 300 count subsamples for each sample type, then combine only for BIOMASS and COUNT. Start with raw data
rawZoop.300 <- rarifyCounts(rawZp.1, sampID = c('UID', 'SAMPLE_TYPE'),
                            abund = 'ABUNDANCE_TOTAL',
                            subsize = 300,
                            seed = 1234)

zpCts.300 <- convertZoopCts_NLA(rawZoop.300, sampID='UID', sampType='SAMPLE_TYPE',
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

zpCts.all <- rbind(zonwCts, zpCts) %>%
  filter(COUNT>0) %>%
  merge(taxa, by = c('TAXA_ID'))

zpCts.all.dist <- assignDistinct(zpCts.all, sampleID = c('UID', 'SAMPLE_TYPE'), taxlevels=c('PHYLUM','CLASS','SUBCLASS','ORDER','SUBORDER','FAMILY','GENUS','SPECIES','SUBSPECIES'))

zpCts.all.300 <- rbind(zonwCts.300, zpCts.300) %>%
  filter(COUNT_300>0) %>%
  merge(taxa, by = c('TAXA_ID'))

zpCts.all.300.dist <- assignDistinct(zpCts.all.300, sampleID = c('UID', 'SAMPLE_TYPE'), taxlevels=c('PHYLUM','CLASS','SUBCLASS','ORDER','SUBORDER','FAMILY','GENUS','SPECIES','SUBSPECIES')) %>%
  plyr::rename(c('IS_DISTINCT' = 'IS_DISTINCT_300'))



