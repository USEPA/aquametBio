# Test_ZooplanktonCode.r
# Use this code to pull data and test functions to make sure the outputs
# match up with values in database.

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

rawZp <- dbGet('ALL_THE_NLA', 'tblZOOPRAW', where = "PARAMETER IN('ABUNDANCE_TOTAL', 'BIOMASS_FACTOR', 'CONCENTRATED_VOLUME', 'VOLUME_COUNTED') AND
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

# Compare to values in database in tblZOOPCNT
curCts <- dbGet('ALL_THE_NLA', 'tblZOOPCNT', where = "PARAMETER IN('BIOMASS', 'COUNT', 'DENSITY') AND SAMPLE_TYPE IN('ZOCN', 'ZOFN')")

compareCts <- pivot_longer(zpCts, cols=COUNT:DENSITY, names_to = 'PARAMETER', values_to='RESULT') %>%
  merge(curCts, by=c('UID', 'SAMPLE_TYPE', 'TAXA_ID', 'PARAMETER'), all=TRUE) %>%
  mutate(RESULT.x=as.numeric(RESULT.x), RESULT.y=as.numeric(RESULT.y))

diffCts <- filter(compareCts, abs(RESULT.x-RESULT.y)>0.01) %>%
  mutate(diffFactor = RESULT.x/RESULT.y) # 0 differences now

msgCur <- filter(compareCts, is.na(RESULT.y))
msgNew <- filter(compareCts, is.na(RESULT.x))

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
msgNew.zonw <- filter(compareZONW, is.na(RESULT.x))



taxa <- dbGet('ALL_THE_NLA', 'tblTAXA', where = "ASSEMBLAGE_NAME='ZOOPLANKTON'") %>%
  pivot_wider(id_cols = c('TAXA_ID'), names_from='PARAMETER', values_from='RESULT')
