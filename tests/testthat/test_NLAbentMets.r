library(testthat)
library(aquametBio)

context("Test NLA benthic metric and MMI calculations")


test_that("Benthic Taxonomy metric values correct",
          {
            testOut <- calcBentTaxMets(bentctsNLA_test,inTaxa=bentTaxa_nla,sampID=c('SITE_ID')
                                        ,dist='IS_DISTINCT',ct='TOTAL')
            varLong <- names(testOut)[names(testOut) %nin% c('SITE_ID')]
            testOut.long <- reshape(testOut, idvar = c('SITE_ID'), direction = 'long',
                                    varying = varLong,timevar = 'PARAMETER',
                                    v.names = 'RESULT', times = varLong)
            # testOut.long <- data.table::melt(testOut,id.vars=c('SITE_ID')
            #                                ,variable.name='PARAMETER',value.name='RESULT') %>%
            #   plyr::mutate(PARAMETER=as.character(PARAMETER))
            compOut <- merge(bentMetsNLA_test,testOut.long,by=c('SITE_ID','PARAMETER'))
            expect_true(nrow(compOut)==550)
            expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.0001)

          })

test_that("Benthic FFG metric values correct",
{
  testOut <- calcBentFFGmets(bentctsNLA_test,inTaxa=bentTaxa_nla,sampID=c('SITE_ID')
                              ,dist='IS_DISTINCT',ct='TOTAL',ffg='FFG')
  varLong <- names(testOut)[names(testOut) %nin% c('SITE_ID')]
  testOut.long <- reshape(testOut, idvar = c('SITE_ID'), direction = 'long',
                          varying = varLong,timevar = 'PARAMETER',
                          v.names = 'RESULT', times = varLong)
  # testOut.long <- data.table::melt(testOut,id.vars=c('SITE_ID')
  #                                ,variable.name='PARAMETER',value.name='RESULT') %>%
  #   plyr::mutate(PARAMETER=as.character(PARAMETER))
  compOut <- merge(bentMetsNLA_test,testOut.long,by=c('SITE_ID','PARAMETER'))
  expect_true(nrow(compOut)==180)
  expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.0001)

})


test_that("Benthic Habit metric values correct",
{
  testOut <- calcBentHabitMets(bentctsNLA_test,inTaxa=bentTaxa_nla,sampID=c('SITE_ID')
                         ,dist='IS_DISTINCT',ct='TOTAL',habit='HABIT')
  varLong <- names(testOut)[names(testOut) %nin% c('SITE_ID')]
  testOut.long <- reshape(testOut, idvar = c('SITE_ID'), direction = 'long',
                          varying = varLong,timevar = 'PARAMETER',
                          v.names = 'RESULT', times = varLong)
  # testOut.long <- data.table::melt(testOut,id.vars=c('SITE_ID')
  #                                ,variable.name='PARAMETER',value.name='RESULT') %>%
  #   plyr::mutate(PARAMETER=as.character(PARAMETER))
  compOut <- merge(bentMetsNLA_test,testOut.long,by=c('SITE_ID','PARAMETER'))
  expect_true(nrow(compOut)==150)
  expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.0001)

})

test_that("Benthic tolerance metric values correct",
{
  testOut <- calcBentTolMets(bentctsNLA_test,inTaxa=bentTaxa_nla,sampID=c('SITE_ID')
                           ,dist='IS_DISTINCT',ct='TOTAL',ptv='PTV')
  varLong <- names(testOut)[names(testOut) %nin% c('SITE_ID')]
  testOut.long <- reshape(testOut, idvar = c('SITE_ID'), direction = 'long',
                          varying = varLong,timevar = 'PARAMETER',
                          v.names = 'RESULT', times = varLong)
  # testOut.long <- data.table::melt(testOut,id.vars=c('SITE_ID')
  #                                ,variable.name='PARAMETER',value.name='RESULT') %>%
  #   plyr::mutate(PARAMETER=as.character(PARAMETER))
  compOut <- merge(bentMetsNLA_test,testOut.long,by=c('SITE_ID','PARAMETER'))
  expect_true(nrow(compOut)==280)
  expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.0001)

})


test_that("Benthic dominance and diversity metric values correct",
{
  testOut <- calcBentDominMets(bentctsNLA_test,inTaxa=bentTaxa_nla,sampID=c('SITE_ID')
                         ,dist='IS_DISTINCT',ct='TOTAL')
  varLong <- names(testOut)[names(testOut) %nin% c('SITE_ID')]
  testOut.long <- reshape(testOut, idvar = c('SITE_ID'), direction = 'long',
                          varying = varLong,timevar = 'PARAMETER',
                          v.names = 'RESULT', times = varLong)
  # testOut.long <- data.table::melt(testOut,id.vars=c('SITE_ID')
  #                                ,variable.name='PARAMETER',value.name='RESULT') %>%
  #   plyr::mutate(PARAMETER=as.character(PARAMETER))
  compOut <- merge(bentMetsNLA_test,testOut.long,by=c('SITE_ID','PARAMETER'))
  expect_true(nrow(compOut)==70)
  expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.0001)
  })

test_that("All benthic metric values correct",
          {
            testOut <- calcAllBentMets(indf=bentctsNLA_test,inTaxa=bentTaxa_nla,sampID=c('SITE_ID')
                                       ,dist='IS_DISTINCT',ct='TOTAL',taxa_id='TAXA_ID',ffg='FFG'
                                       ,habit='HABIT',ptv='PTV')
            varLong <- names(testOut)[names(testOut) %nin% c('SITE_ID')]
            testOut.long <- reshape(testOut, idvar = c('SITE_ID'), direction = 'long',
                                    varying = varLong,timevar = 'PARAMETER',
                                    v.names = 'RESULT', times = varLong)
            # testOut.long <- data.table::melt(testOut,id.vars=c('SITE_ID')
            #                                ,variable.name='PARAMETER',value.name='RESULT') %>%
            #   plyr::mutate(PARAMETER=as.character(PARAMETER))
            compOut <- merge(bentMetsNLA_test,testOut.long,by=c('SITE_ID','PARAMETER'))
            expect_true(nrow(compOut)==1250)
            expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.0001)
          })

indf.eco <- subset(bentMetsNLA_test,select=c('SITE_ID','ECO_BIO')) %>%
  unique() %>% merge(bentctsNLA_test,by='SITE_ID')

test_that("MMI metrics correct",
          {
            testOut <- calcNLA_BentMMImets(inCts=indf.eco,inTaxa=bentTaxa_nla,sampID=c('SITE_ID')
                                       ,dist='IS_DISTINCT',ct='TOTAL',taxa_id='TAXA_ID',ffg='FFG'
                                       ,habit='HABIT',ptv='PTV',ecoreg='ECO_BIO')
            varLong <- names(testOut)[names(testOut) %nin% c('SITE_ID','ECO_BIO')]
            testOut.long <- reshape(testOut, idvar = c('SITE_ID'), direction = 'long',
                                    varying = varLong,timevar = 'PARAMETER',
                                    v.names = 'RESULT', times = varLong)
            # testOut.long <- data.table::melt(testOut,id.vars=c('SITE_ID','ECO_BIO')
            #                                ,variable.name='PARAMETER',value.name='RESULT',na.rm=T) %>%
            #   plyr::mutate(PARAMETER=as.character(PARAMETER))
            testOut.long <- subset(testOut.long, !is.na(RESULT))
            compOut <- merge(bentMetsNLA_test,testOut.long,by=c('SITE_ID','PARAMETER'))
            expect_true(nrow(compOut)==70)
            expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.0001)

          })

testIn <- reshape(bentMetsNLA_test, idvar = c('SITE_ID','ECO_BIO'), direction = 'wide',
                 timevar = 'PARAMETER', v.names = 'RESULT')
names(testIn) <- gsub('RESULT\\.', '', names(testIn))
# testIn <- data.table::dcast(bentMetsNLA_test,SITE_ID+ECO_BIO~PARAMETER,value.var='RESULT')

test_that("NRSA Benthic MMI scores correct",
{
  testOut <- calcNLA_BenthicMMI(testIn,sampID=c('SITE_ID')
                                 ,ecoreg='ECO_BIO',totlnind='TOTLNIND')
  varLong <- names(testOut)[names(testOut) %nin% c('SITE_ID','ECO_BIO')]
  testOut.long <- reshape(testOut, idvar = c('SITE_ID'), direction = 'long',
                          varying = varLong,timevar = 'PARAMETER',
                          v.names = 'RESULT', times = varLong)
  testOut.long <- subset(testOut.long, !is.na(RESULT))

  varLong.1 <- names(bentMMI_NLA_test)[names(bentMMI_NLA_test) %nin% c('SITE_ID','ECO_BIO','TOTLNIND')]
  bentMMI_NLA_test.long <- reshape(bentMMI_NLA_test, idvar = 'SITE_ID', direction = 'long',
                                   varying = varLong.1, times = varLong.1, timevar = 'PARAMETER',
                                   v.names = 'RESULT')
  bentMMI_NLA_test.long <- subset(bentMMI_NLA_test.long, !is.na(RESULT))
  print(unique(bentMMI_NLA_test.long$PARAMETER))
  print(nrow(bentMMI_NLA_test.long))
  print(unique(testOut.long$PARAMETER))
  print(nrow(testOut.long))
  compOut <- merge(bentMMI_NLA_test.long, testOut.long, by=c('SITE_ID','PARAMETER'))
  # expect_true(nrow(compOut)==73)
  compOut.cond <- subset(compOut, PARAMETER=='BENT_MMI_COND')
  expect_equal(compOut.cond$RESULT.x,compOut.cond$RESULT.y)
  compOut.res <- subset(compOut,PARAMETER!='BENT_MMI_COND')
  compOut.res$RESULT.x <- as.numeric(compOut.res$RESULT.x)
  compOut.res$RESULT.y <- as.numeric(compOut.res$RESULT.y)
  expect_equal(compOut.res$RESULT.x,compOut.res$RESULT.y,tolerance=0.01)
  expect_equal(nrow(compOut.res), 63)

})

