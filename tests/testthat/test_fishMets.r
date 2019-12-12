library(testthat)
library(aquametBio)

context("Test fish metric and MMI calculations")

test_that("Fish tolerance metric values correct",
{
  testOut <- calcFishTolMets(fishCts_test,sampID=c('UID'),dist='IS_DISTINCT',ct='FINAL_CT'
                             ,taxa_id='TAXA_ID',tol='TOLERANCE_NRSA',tolval='TOL_VAL_EMAPW'
                             ,vel='VEL_NRSA',habitat='HABITAT_NRSA',trophic='TROPHIC_NRSA'
                             ,migr='MIGR_NRSA',nonnat='NON_NATIVE')
  varLong <- names(testOut)[names(testOut) %nin% c('UID')]
  testOut.long <- reshape(testOut, idvar = c('UID'), direction = 'long',
                          varying = varLong,timevar = 'PARAMETER',
                          v.names = 'RESULT', times = varLong)
  testOut.long$RESULT <- as.numeric(testOut.long$RESULT)
  # testOut.long <- data.table::melt(testOut,id.vars=c('UID')
  #                                ,variable.name='PARAMETER',value.name='RESULT') %>%
  #   plyr::mutate(PARAMETER=as.character(PARAMETER),RESULT=as.numeric(RESULT))
  compOut <- merge(fishMet_test,testOut.long,by=c('UID','PARAMETER'))
  expect_true(nrow(compOut)==740)
  expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.0001)

})

# Run again with missing traits
test_that("Fish tolerance metric values without native or habitat correct",
{
  testOut <- calcFishTolMets(fishCts_test,sampID=c('UID'),dist='IS_DISTINCT',ct='FINAL_CT'
                             ,taxa_id='TAXA_ID',tol='TOLERANCE_NRSA',tolval='TOL_VAL_EMAPW'
                             ,vel='VEL_NRSA',trophic='TROPHIC_NRSA'
                             ,migr='MIGR_NRSA')
  varLong <- names(testOut)[names(testOut) %nin% c('UID')]
  testOut.long <- reshape(testOut, idvar = c('UID'), direction = 'long',
                          varying = varLong,timevar = 'PARAMETER',
                          v.names = 'RESULT', times = varLong)
  testOut.long$RESULT <- as.numeric(testOut.long$RESULT)
  # testOut.long <- data.table::melt(testOut,id.vars=c('UID')
  #                                ,variable.name='PARAMETER',value.name='RESULT') %>%
  #   plyr::mutate(PARAMETER=as.character(PARAMETER),RESULT=as.numeric(RESULT))
  compOut <- merge(fishMet_test,testOut.long,by=c('UID','PARAMETER'))
  expect_true(nrow(compOut)==330)
  expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.0001)

})

test_that("Fish taxonomic metric values correct",
{
  testOut <- calcFishTaxMets(fishCts_test,fishTaxa,sampID=c('UID'),dist='IS_DISTINCT',ct='FINAL_CT'
                             ,taxa_id='TAXA_ID',nonnat='NON_NATIVE',family='FAMILY'
                             ,genus='GENUS',comname='FINAL_NAME')
  varLong <- names(testOut)[names(testOut) %nin% c('UID')]
  testOut.long <- reshape(testOut, idvar = c('UID'), direction = 'long',
                          varying = varLong,timevar = 'PARAMETER',
                          v.names = 'RESULT', times = varLong)
  testOut.long$RESULT <- as.numeric(testOut.long$RESULT)
  # testOut.long <- data.table::melt(testOut,id.vars=c('UID')
  #                                ,variable.name='PARAMETER',value.name='RESULT') %>%
  #   plyr::mutate(PARAMETER=as.character(PARAMETER),RESULT=as.numeric(RESULT))
  compOut <- merge(fishMet_test,testOut.long,by=c('UID','PARAMETER'))
  expect_true(nrow(compOut)==420)
  expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.0001)

})

test_that("Fish trophic metric values correct",
{
  testOut <- calcFishTrophicMets(fishCts_test,fishTaxa,sampID=c('UID'),dist='IS_DISTINCT',ct='FINAL_CT'
                             ,taxa_id='TAXA_ID',nonnat='NON_NATIVE',trophic='TROPHIC_NRSA'
                             ,habitat='HABITAT_NRSA')
  varLong <- names(testOut)[names(testOut) %nin% c('UID')]
  testOut.long <- reshape(testOut, idvar = c('UID'), direction = 'long',
                          varying = varLong,timevar = 'PARAMETER',
                          v.names = 'RESULT', times = varLong)
  testOut.long$RESULT <- as.numeric(testOut.long$RESULT)
  # testOut.long <- data.table::melt(testOut,id.vars=c('UID')
  #                                ,variable.name='PARAMETER',value.name='RESULT') %>%
  #   plyr::mutate(PARAMETER=as.character(PARAMETER),RESULT=as.numeric(RESULT))
  compOut <- merge(fishMet_test,testOut.long,by=c('UID','PARAMETER'))
  expect_true(nrow(compOut)==360)
  expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.0001)

})

test_that("Fish other metric values correct",
{
  testOut <- calcFishOtherMets(fishCts_test,fishTaxa,sampID=c('UID'),dist='IS_DISTINCT',ct='FINAL_CT'
                               ,taxa_id='TAXA_ID',vel='VEL_NRSA'
                               ,migr='MIGR_NRSA', reprod='REPROD_NRSA', temp='TEMP_NRSA'
                               ,nonnat='NON_NATIVE')
  varLong <- names(testOut)[names(testOut) %nin% c('UID')]
  testOut.long <- reshape(testOut, idvar = c('UID'), direction = 'long',
                          varying = varLong,timevar = 'PARAMETER',
                          v.names = 'RESULT', times = varLong)
  testOut.long$RESULT <- as.numeric(testOut.long$RESULT)
  # testOut.long <- data.table::melt(testOut,id.vars=c('UID')
  #                                ,variable.name='PARAMETER',value.name='RESULT') %>%
  #   plyr::mutate(PARAMETER=as.character(PARAMETER),RESULT=as.numeric(RESULT))
  compOut <- merge(fishMet_test,testOut.long,by=c('UID','PARAMETER'),all.y=T)
  expect_true(nrow(compOut)==360)
  expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.0001)

})

test_that("Fish native/alien metric values correct",
{
  testOut <- calcFishNativeMets(fishCts_test,sampID=c('UID'),dist='IS_DISTINCT',ct='FINAL_CT'
                               ,taxa_id='TAXA_ID',nonnat='NON_NATIVE')
  varLong <- names(testOut)[names(testOut) %nin% c('UID')]
  testOut.long <- reshape(testOut, idvar = c('UID'), direction = 'long',
                          varying = varLong,timevar = 'PARAMETER',
                          v.names = 'RESULT', times = varLong)
  testOut.long$RESULT <- as.numeric(testOut.long$RESULT)
  # testOut.long <- data.table::melt(testOut,id.vars=c('UID')
  #                                ,variable.name='PARAMETER',value.name='RESULT') %>%
  #   plyr::mutate(PARAMETER=as.character(PARAMETER),RESULT=as.numeric(RESULT))
  compOut <- merge(fishMet_test,testOut.long,by=c('UID','PARAMETER'),all.y=T)
  expect_true(nrow(compOut)==90)
  expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.0001)

})

test_that("Fish anomaly metric values correct",
{
  testOut <- calcFishAnomMets(fishCts_test,sampID=c('UID'),ct='FINAL_CT'
                                ,anomct='ANOM_CT')
  varLong <- names(testOut)[names(testOut) %nin% c('UID')]
  testOut.long <- reshape(testOut, idvar = c('UID'), direction = 'long',
                          varying = varLong,timevar = 'PARAMETER',
                          v.names = 'RESULT', times = varLong)
  testOut.long$RESULT <- as.numeric(testOut.long$RESULT)
  # testOut.long <- data.table::melt(testOut,id.vars=c('UID')
  #                                ,variable.name='PARAMETER',value.name='RESULT') %>%
  #   plyr::mutate(PARAMETER=as.character(PARAMETER),RESULT=as.numeric(RESULT))
  compOut <- merge(fishMet_test,testOut.long,by=c('UID','PARAMETER'))
  expect_true(nrow(compOut)==10)
  expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.0001)

})

test_that("All fish metric values correct",
{
  testOut <- calcAllFishMets(fishCts_test,fishTaxa,sampID=c('UID'),dist='IS_DISTINCT',
                             ct='FINAL_CT',anomct='ANOM_CT',taxa_id='TAXA_ID',
                             tol='TOLERANCE_NRSA',tolval='TOL_VAL_EMAPW',vel='VEL_NRSA',
                             habitat='HABITAT_NRSA',trophic='TROPHIC_NRSA',
                             migr='MIGR_NRSA',nonnat='NON_NATIVE',
                             reprod='REPROD_NRSA', temp='TEMP_NRSA',
                             family='FAM_OR_CLS',genus='GENUS',comname='FINAL_NAME')
  varLong <- names(testOut)[names(testOut) %nin% c('UID')]
  testOut.long <- reshape(testOut, idvar = c('UID'), direction = 'long',
                          varying = varLong,timevar = 'PARAMETER',
                          v.names = 'RESULT', times = varLong)
  # testOut.long <- data.table::melt(testOut,id.vars=c('UID')
  #                                ,variable.name='PARAMETER',value.name='RESULT') %>%
  #   plyr::mutate(PARAMETER=as.character(PARAMETER),RESULT=as.numeric(RESULT))
  compOut <- merge(fishMet_test,testOut.long,by=c('UID','PARAMETER'))
  expect_true(nrow(compOut)==1740)
  expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.0001)
})


ecoTest <- data.frame(UID=c(11222, 11703, 11711, 12384, 14080, 14275, 14315, 15194, 15196, 15198)
                      ,AGGR_ECO9_2015=c('UMW','CPL','CPL','CPL','CPL','SAP','TPL','XER','XER','WMT'),stringsAsFactors=F)
fishCts.eco <- merge(fishCts_test,ecoTest,by='UID')

test_that("MMI fish metric values correct",
{
  testOut <- calcNRSA_FishMMImets(fishCts.eco,fishTaxa,sampID=c('UID'),dist='IS_DISTINCT',
                             ct='FINAL_CT',taxa_id='TAXA_ID',
                             tol='TOLERANCE_NRSA',vel='VEL_NRSA',
                             habitat='HABITAT_NRSA',trophic='TROPHIC_NRSA',
                             migr='MIGR_NRSA',nonnat='NON_NATIVE',
                             reprod='REPROD_NRSA',
                             family='FAMILY',genus='GENUS',comname='FINAL_NAME',
                             ecoreg='AGGR_ECO9_2015')
  varLong <- names(testOut)[names(testOut) %nin% c('UID','AGGR_ECO9_2015')]
  testOut.long <- reshape(testOut, idvar = c('UID','AGGR_ECO9_2015'), direction = 'long',
                          varying = varLong,timevar = 'PARAMETER',
                          v.names = 'RESULT', times = varLong)
  testOut.long <- subset(testOut.long, !is.na(RESULT))
  testOut.long$RESULT <- as.numeric(testOut.long$RESULT)
  # testOut.long <- data.table::melt(testOut,id.vars=c('UID','AGGR_ECO9_2015')
  #                                ,variable.name='PARAMETER',value.name='RESULT',na.rm=T) %>%
  #   plyr::mutate(PARAMETER=as.character(PARAMETER),RESULT=as.numeric(RESULT))
  compOut <- merge(fishMet_test,testOut.long,by=c('UID','PARAMETER'))
  expect_true(nrow(compOut)==90)
  expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.0001)
})

# Create data frame of log10(watershed area) values

lwsarea <- data.frame(UID=c(11222, 11703, 11711, 12384, 14080, 14275, 14315, 15194, 15196, 15198)
                      ,LWSAREA=c(3.88288,6.36191,6.36425,3.81217,-1.02872,4.41280,2.80386,2.28921,2.45688,2.02662))
fishws_test <- merge(fishMet_test,lwsarea,by='UID') %>% merge(ecoTest,by='UID')
# Put fish metrics in wide format
fishws_test.1 <- reshape(fishws_test, idvar = c('UID','SAMPLE_TYPE','LWSAREA','AGGR_ECO9_2015'), direction = 'wide',
                         v.names = 'RESULT', timevar = 'PARAMETER')
names(fishws_test.1) <- gsub("RESULT\\.", "", names(fishws_test.1))
# fishws_test.1 <- data.table::dcast(fishws_test,UID+SAMPLE_TYPE+LWSAREA+AGGR_ECO9_2015~PARAMETER,value.var='RESULT')

test_that("Fish MMI scores correct",
          {
            testOut <- calcFishMMI(fishws_test.1,sampID=c('UID','SAMPLE_TYPE'),ecoreg='AGGR_ECO9_2015'
                                   ,lwsarea='LWSAREA')
            testOut.long <- reshape(testOut, idvar = c('UID','AGGR_ECO9_2015'), direction = 'long',
                                    varying = 'MMI_FISH', timevar = 'PARAMETER',
                                    v.names = 'RESULT', times = 'MMI_FISH')
            testOut.long <- subset(testOut.long, !is.na(RESULT))
            # testOut.long <- data.table::melt(testOut,id.vars=c('UID','AGGR_ECO9_2015','SAMPLE_TYPE')
            #                                ,variable.name='PARAMETER',value.name='RESULT',na.rm=T) %>%
            #   plyr::mutate(PARAMETER=as.character(PARAMETER))
            fishMMI_test.long <- reshape(fishMMI_test, idvar = c('UID'), direction = 'long',
                                         varying = 'MMI_FISH', times = 'MMI_FISH', timevar = 'PARAMETER',
                                         v.names = 'RESULT')
            # fishMMI_test.long <- data.table::melt(fishMMI_test,id.vars=c('UID'),measure.vars='MMI_FISH'
            #                                     ,variable.name='PARAMETER',value.name='RESULT')
            compOut <- merge(fishMMI_test.long,testOut.long,by=c('UID','PARAMETER'))
            expect_true(nrow(compOut)==10)
            expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.0001)

          })

# Create wsarea from lwsarea above for samples in fishws_test.1
fishws_test.2 <- fishws_test.1
fishws_test.2$WSAREA <- with(fishws_test.2, 10^LWSAREA)
fishws_test.2 <- merge(fishws_test.1, fishMMI_test[,c('UID','MMI_FISH')], by = 'UID')
# fishws_test.2 <- plyr::mutate(fishws_test.1,WSAREA=10^LWSAREA) %>%
#   merge(fishMMI_test[,c('UID','MMI_FISH')],by='UID')

test_that("Fish MMI scores correct",
          {
            testOut <- assignFishCondition(fishws_test.2,sampID=c('UID','SAMPLE_TYPE'),ecoreg='AGGR_ECO9_2015'
                                   ,wsarea='WSAREA',totlnind='TOTLNIND',mmi='MMI_FISH')
            testOut.1 <- testOut[,c('UID','SAMPLE_TYPE','MMI_FISH','FISH_MMI_COND')]
            testOut.long <- reshape(testOut.1, idvar = c('UID', 'SAMPLE_TYPE'), direction = 'long',
                                    varying = c('MMI_FISH','FISH_MMI_COND'), times = c('MMI_FISH','FISH_MMI_COND'),
                                    timevar = 'PARAMETER', v.names = 'RESULT')
            # testOut.long <- data.table::melt(testOut,id.vars=c('UID','SAMPLE_TYPE')
            #                                ,variable.name='PARAMETER',value.name='RESULT',na.rm=T) %>%
            #   plyr::mutate(PARAMETER=as.character(PARAMETER))
            fishMMI_test.long <- reshape(fishMMI_test, idvar = 'UID', direction = 'long',
                                         varying = 'FISH_MMI_COND', times = 'FISH_MMI_COND',
                                         v.names = 'RESULT', timevar = 'PARAMETER')
            # fishMMI_test.long <- data.table::melt(fishMMI_test,id.vars=c('UID'),measure.vars='FISH_MMI_COND'
            #                                     ,variable.name='PARAMETER',value.name='RESULT')
            compOut <- merge(fishMMI_test.long,testOut.long,by=c('UID','PARAMETER'))
            expect_true(nrow(compOut)==10)
            expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.0001)

          })
