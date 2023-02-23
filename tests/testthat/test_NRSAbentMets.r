# context("Test NRSA benthic metric and MMI calculations")


test_that("Data prep correct",
{
  testOut <- prepBentCts_WSA(inCts=bentCts_test,inTaxa=bentTaxa_nrsa
                             ,sampID=c('UID','SAMPLE_TYPE','SAMPLE_CAT')
                             ,ct='TOTAL300'
                             ,taxa_id='TAXA_ID')
  expect_equal(nrow(testOut),nrow(indf_test))
  expect_equal(names(testOut),names(indf_test))
  testOut.long <- reshape(testOut, idvar = c('UID','SAMPLE_TYPE','SAMPLE_CAT','TAXA_ID'), direction = 'long',
                          varying = c('TOTAL','IS_DISTINCT'),timevar = 'variable', v.names = 'value',
                          times = c('TOTAL','IS_DISTINCT'))
#  testOut.long <- data.table::melt(testOut,id.vars=c('UID','SAMPLE_TYPE','SAMPLE_CAT','TAXA_ID'))
  indf.long <- reshape(indf_test, idvar = c('UID','SAMPLE_TYPE','SAMPLE_CAT','TAXA_ID'), direction = 'long',
                       varying = c('TOTAL','IS_DISTINCT'),timevar = 'variable', v.names = 'value',
                       times = c('TOTAL','IS_DISTINCT'))
#  indf.long <- data.table::melt(indf_test,id.vars=c('UID','SAMPLE_TYPE','SAMPLE_CAT','TAXA_ID'))
  compOut <- merge(testOut.long,indf.long,by=c('UID','SAMPLE_TYPE','SAMPLE_CAT','TAXA_ID','variable'))
  expect_equal(compOut$value.x,compOut$value.y)
})


test_that("Benthic Taxonomy metric values correct",
          {
            testOut <- calcBentTaxMets(indf_test,inTaxa=bentTaxa_nrsa,sampID=c('UID','SAMPLE_TYPE','SAMPLE_CAT')
                                        ,dist='IS_DISTINCT',ct='TOTAL')
            varLong <- names(testOut)[names(testOut) %nin% c('UID','SAMPLE_TYPE','SAMPLE_CAT')]
            testOut.long <- reshape(testOut, idvar = c('UID','SAMPLE_TYPE','SAMPLE_CAT'), direction = 'long',
                                    varying = varLong,timevar = 'PARAMETER',
                                    v.names = 'RESULT', times = varLong)
            # testOut.long <- data.table::melt(testOut,id.vars=c('UID','SAMPLE_TYPE','SAMPLE_CAT')
            #                                ,variable.name='PARAMETER',value.name='RESULT') %>%
            #   plyr::mutate(PARAMETER=as.character(PARAMETER))
            compOut <- merge(bentMet_test,testOut.long,by=c('UID','SAMPLE_TYPE','SAMPLE_CAT','PARAMETER'))
            expect_true(nrow(compOut)==550)
            expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.0001)

          })

test_that("Benthic FFG metric values correct",
{
  testOut <- calcBentFFGmets(indf_test,inTaxa=bentTaxa_nrsa,sampID=c('UID','SAMPLE_TYPE','SAMPLE_CAT')
                              ,dist='IS_DISTINCT',ct='TOTAL',ffg='FFG_WSA')
  varLong <- names(testOut)[names(testOut) %nin% c('UID','SAMPLE_TYPE','SAMPLE_CAT')]
  testOut.long <- reshape(testOut, idvar = c('UID','SAMPLE_TYPE','SAMPLE_CAT'), direction = 'long',
                          varying = varLong,timevar = 'PARAMETER',
                          v.names = 'RESULT', times = varLong)
  # testOut.long <- data.table::melt(testOut,id.vars=c('UID','SAMPLE_TYPE','SAMPLE_CAT')
  #                                ,variable.name='PARAMETER',value.name='RESULT') %>%
  #   plyr::mutate(PARAMETER=as.character(PARAMETER))
  compOut <- merge(bentMet_test,testOut.long,by=c('UID','SAMPLE_TYPE','SAMPLE_CAT','PARAMETER'))
  expect_true(nrow(compOut)==180)
  expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.0001)

})


test_that("Benthic Habit metric values correct",
{
  testOut <- calcBentHabitMets(indf_test,inTaxa=bentTaxa_nrsa,sampID=c('UID','SAMPLE_TYPE','SAMPLE_CAT')
                         ,dist='IS_DISTINCT',ct='TOTAL',habit='HABIT_WSA')
  varLong <- names(testOut)[names(testOut) %nin% c('UID','SAMPLE_TYPE','SAMPLE_CAT')]
  testOut.long <- reshape(testOut, idvar = c('UID','SAMPLE_TYPE','SAMPLE_CAT'), direction = 'long',
                          varying = varLong,timevar = 'PARAMETER',
                          v.names = 'RESULT', times = varLong)

  # testOut.long <- data.table::melt(testOut,id.vars=c('UID','SAMPLE_TYPE','SAMPLE_CAT')
  #                                ,variable.name='PARAMETER',value.name='RESULT') %>%
  #   plyr::mutate(PARAMETER=as.character(PARAMETER))
  compOut <- merge(bentMet_test,testOut.long,by=c('UID','SAMPLE_TYPE','SAMPLE_CAT','PARAMETER'))
  expect_true(nrow(compOut)==150)
  expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.0001)

})

test_that("Benthic tolerance metric values correct",
{
  testOut <- calcBentTolMets(indf_test,inTaxa=bentTaxa_nrsa,sampID=c('UID','SAMPLE_TYPE','SAMPLE_CAT')
                           ,dist='IS_DISTINCT',ct='TOTAL',ptv='PTV_WSA')
  varLong <- names(testOut)[names(testOut) %nin% c('UID','SAMPLE_TYPE','SAMPLE_CAT')]
  testOut.long <- reshape(testOut, idvar = c('UID','SAMPLE_TYPE','SAMPLE_CAT'), direction = 'long',
                          varying = varLong,timevar = 'PARAMETER',
                          v.names = 'RESULT', times = varLong)
  # testOut.long <- data.table::melt(testOut,id.vars=c('UID','SAMPLE_TYPE','SAMPLE_CAT')
  #                                ,variable.name='PARAMETER',value.name='RESULT') %>%
  #   plyr::mutate(PARAMETER=as.character(PARAMETER))
  compOut <- merge(bentMet_test,testOut.long,by=c('UID','SAMPLE_TYPE','SAMPLE_CAT','PARAMETER'))
  expect_true(nrow(compOut)==280)
  expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.0001)

})


test_that("Benthic dominance and diversity metric values correct",
{
  testOut <- calcBentDominMets(indf_test,inTaxa=bentTaxa_nrsa,sampID=c('UID','SAMPLE_TYPE','SAMPLE_CAT')
                         ,dist='IS_DISTINCT',ct='TOTAL')
  varLong <- names(testOut)[names(testOut) %nin% c('UID','SAMPLE_TYPE','SAMPLE_CAT')]
  testOut.long <- reshape(testOut, idvar = c('UID','SAMPLE_TYPE','SAMPLE_CAT'), direction = 'long',
                          varying = varLong,timevar = 'PARAMETER',
                          v.names = 'RESULT', times = varLong)
  # testOut.long <- data.table::melt(testOut,id.vars=c('UID','SAMPLE_TYPE','SAMPLE_CAT')
  #                                ,variable.name='PARAMETER',value.name='RESULT') %>%
  #   plyr::mutate(PARAMETER=as.character(PARAMETER))
  compOut <- merge(bentMet_test,testOut.long,by=c('UID','SAMPLE_TYPE','SAMPLE_CAT','PARAMETER'))
  expect_true(nrow(compOut)==70)
  expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.0001)
  })

test_that("All benthic metric values correct",
          {
            testOut <- calcAllBentMets(indf=indf_test,inTaxa=bentTaxa_nrsa,sampID=c('UID','SAMPLE_TYPE','SAMPLE_CAT')
                                       ,dist='IS_DISTINCT',ct='TOTAL',taxa_id='TAXA_ID',ffg='FFG_WSA'
                                       ,habit='HABIT_WSA',ptv='PTV_WSA')
            varLong <- names(testOut)[names(testOut) %nin% c('UID','SAMPLE_TYPE','SAMPLE_CAT')]
            testOut.long <- reshape(testOut, idvar = c('UID','SAMPLE_TYPE','SAMPLE_CAT'), direction = 'long',
                                    varying = varLong,timevar = 'PARAMETER',
                                    v.names = 'RESULT', times = varLong)
            # testOut.long <- data.table::melt(testOut,id.vars=c('UID','SAMPLE_TYPE','SAMPLE_CAT')
            #                                ,variable.name='PARAMETER',value.name='RESULT') %>%
            #   plyr::mutate(PARAMETER=as.character(PARAMETER))
            compOut <- merge(bentMet_test,testOut.long,by=c('UID','SAMPLE_TYPE','SAMPLE_CAT','PARAMETER'))
            expect_true(nrow(compOut)==1250)
            expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.0001)
          })



test_that("MMI metrics correct",
          {
            ecoTest <- data.frame(UID=c(11245,11703,11821,12483,12571,12999,13971,14822,15073,15803)
                                  ,AGGR_ECO9_2015=c('TPL','CPL','NAP','SAP','SAP','UMW','CPL','NPL'
                                                    ,'UMW','XER'),stringsAsFactors=F)
            indf.eco <- merge(indf_test,ecoTest,by='UID')

            testOut <- calcNRSA_BentMMImets(inCts=indf.eco,sampID=c('UID','SAMPLE_TYPE','SAMPLE_CAT')
                                       ,dist='IS_DISTINCT',ct='TOTAL',taxa_id='TAXA_ID',ffg='FFG_WSA'
                                       ,habit='HABIT_WSA',ptv='PTV_WSA',ecoreg='AGGR_ECO9_2015')
            varLong <- names(testOut)[names(testOut) %nin% c('UID','SAMPLE_TYPE','SAMPLE_CAT','AGGR_ECO9_2015')]
            testOut.long <- reshape(testOut, idvar = c('UID','SAMPLE_TYPE','SAMPLE_CAT','AGGR_ECO9_2015'), direction = 'long',
                                    varying = varLong, timevar = 'PARAMETER',
                                    v.names = 'RESULT', times = varLong)
            testOut.long <- subset(testOut.long, !is.na(RESULT))
            # testOut.long <- data.table::melt(testOut,id.vars=c('UID','SAMPLE_TYPE','SAMPLE_CAT','AGGR_ECO9_2015')
            #                                ,variable.name='PARAMETER',value.name='RESULT',na.rm=T) %>%
            #   plyr::mutate(PARAMETER=as.character(PARAMETER))
            compOut <- merge(bentMet_test,testOut.long,by=c('UID','SAMPLE_TYPE','SAMPLE_CAT','PARAMETER'))
            expect_true(nrow(compOut)==70)
            expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.0001)

          })


# testIn <- merge(bentMet_test,ecoTest,by='UID') %>%
#   data.table::dcast(UID+SAMPLE_TYPE+SAMPLE_CAT+AGGR_ECO9_2015~PARAMETER,value.var='RESULT')

test_that("NRSA Benthic MMI scores correct",
{
  ecoTest <- data.frame(UID=c(11245,11703,11821,12483,12571,12999,13971,14822,15073,15803)
                        ,AGGR_ECO9_2015=c('TPL','CPL','NAP','SAP','SAP','UMW','CPL','NPL'
                                          ,'UMW','XER'),stringsAsFactors=F)

  testIn <- merge(bentMet_test, ecoTest, by='UID')
  testIn.wide <- reshape(testIn, idvar = c('UID','SAMPLE_TYPE','SAMPLE_CAT','AGGR_ECO9_2015'), direction = 'wide',
                         v.names = 'RESULT', timevar = 'PARAMETER')
  names(testIn.wide) <- gsub('RESULT\\.', '', names(testIn.wide))

  testOut <- calcNRSA_BenthicMMI(testIn.wide,sampID=c('UID','SAMPLE_TYPE','SAMPLE_CAT')
                                 ,ecoreg='AGGR_ECO9_2015',totlnind='TOTLNIND')
  varLong <- names(testOut)[names(testOut) %nin% c('UID','SAMPLE_TYPE','SAMPLE_CAT')]
  testOut.long <- reshape(testOut, idvar = c('UID','AGGR_ECO9_2015','SAMPLE_TYPE','SAMPLE_CAT'), direction = 'long',
                          varying = varLong, timevar = 'PARAMETER',
                          v.names = 'RESULT', times = varLong)
  # testOut.long <- data.table::melt(testOut,id.vars=c('UID','AGGR_ECO9_2015','SAMPLE_TYPE','SAMPLE_CAT')
  #                                ,variable.name='PARAMETER',value.name='RESULT',na.rm=T) %>%
  #   plyr::mutate(PARAMETER=as.character(PARAMETER))
  varLong <- names(bentMMI_test)[names(bentMMI_test) %nin% c('UID','AGGR_ECO9_2015','SAMPLE_TYPE')]
  bentMMI_test.long <- reshape(bentMMI_test, idvar = c('UID'), direction = 'long',
                          varying = varLong, timevar = 'PARAMETER',
                          v.names = 'RESULT', times = varLong)
  # bentMMI_test.long <- data.table::melt(bentMMI_test,id.vars=c('UID')
  #                                     ,variable.name='PARAMETER',value.name='RESULT')
  compOut <- merge(bentMMI_test.long,testOut.long,by=c('UID','PARAMETER'))
  expect_true(nrow(compOut)==80)
  expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.0001)

})
