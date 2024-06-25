# context("Test NLA zooplankton metric and MMI calculations")

test_that("Data prep correct",
          {
            testOut <- convertZoopCts_NLA(rawZoop_test, sampID='UID',
                                          sampType='SAMPLE_TYPE',
                                          rawCt = 'ABUNDANCE_TOTAL',
                                          biofactor = 'BIOMASS_FACTOR',
                                          tow_vol = 'TOW_VOLUME',
                                          vol_ctd = 'VOLUME_COUNTED',
                                          conc_vol = 'CONCENTRATED_VOLUME',
                                          lr_taxon = 'LARGE_RARE_TAXA',
                                          taxa_id = 'TAXA_ID',
                                          subsample=FALSE)
            expect_equal(nrow(testOut), nrow(zpCts_test))
            expect_equal(names(testOut), names(zpCts_test))
            testOut.long <- reshape(testOut, idvar = c('UID','SAMPLE_TYPE','TAXA_ID'), direction = 'long',
                                    varying = c('COUNT', 'BIOMASS', 'DENSITY'),timevar = 'variable', v.names = 'value',
                                    times = c('COUNT', 'BIOMASS', 'DENSITY'))
            cts.long <- reshape(zpCts_test, idvar = c('UID','SAMPLE_TYPE','TAXA_ID'), direction = 'long',
                                 varying = c('COUNT', 'BIOMASS', 'DENSITY'),timevar = 'variable', v.names = 'value',
                                 times = c('COUNT', 'BIOMASS', 'DENSITY'))
            #  indf.long <- data.table::melt(indf_test,id.vars=c('UID','SAMPLE_TYPE','SAMPLE_CAT','TAXA_ID'))
            compOut <- merge(testOut.long, cts.long,by=c('UID','SAMPLE_TYPE','TAXA_ID','variable'))
            expect_equal(compOut$value.x, compOut$value.y)
          })

test_that("Combine fine and coarse mesh samples correct",
          {
           testOut <- prepZoopCombCts_NLA(zpCts_test, sampID = 'UID',
                                          sampType = 'SAMPLE_TYPE',
                                          typeFine = 'ZOFN', typeCoarse = 'ZOCN',
                                          ct = 'COUNT', biomass = 'BIOMASS',
                                          density = 'DENSITY', taxa_id = 'TAXA_ID',
                                          lr_taxon = 'LARGE_RARE_TAXA')
           expect_equal(nrow(testOut), nrow(zonwCts))
           expect_equal(names(testOut), names(zonwCts))
           testOut.long <- reshape(testOut, idvar = c('UID','SAMPLE_TYPE','TAXA_ID'), direction = 'long',
                                   varying = c('COUNT', 'BIOMASS', 'DENSITY'),timevar = 'variable', v.names = 'value',
                                   times = c('COUNT', 'BIOMASS', 'DENSITY'))
           cts.long <- reshape(zpCts_test, idvar = c('UID','SAMPLE_TYPE','TAXA_ID'), direction = 'long',
                               varying = c('COUNT', 'BIOMASS', 'DENSITY'),timevar = 'variable', v.names = 'value',
                               times = c('COUNT', 'BIOMASS', 'DENSITY'))
           compOut <- merge(testOut.long, cts.long,by=c('UID','SAMPLE_TYPE','TAXA_ID','variable'))
           expect_equal(compOut$value.x, compOut$value.y)
          })

test_that("All metrics code correct",
          {
            testOut <- calcZoopAllMets(indata = zonwTest,
                                           inCoarse = zocnTest,
                                           inFine = zofnTest,
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
            compOut <- merge(zoopMets_test, testOut, by=c('UID','PARAMETER'))
            expect_true(nrow(testOut)==nrow(zoopMets_test))
            expect_equal(compOut$RESULT.x, compOut$RESULT.y, tolerance=0.0001)
          })

test_that("MMI scores and condition correct",
          {
            ecoreg <- data.frame(UID = c(2010916, 2010606, 2010513, 2011008, 2010600,
                                         2011314, 2010571, 2011195, 2010632, 2011304),
                                 ECO_BIO = c('EHIGH', 'EHIGH', 'EHIGH', 'UMW', 'EHIGH',
                                             'WMTNS', 'EHIGH', 'PLAINS', 'PLAINS', 'PLAINS'))

            mmiIn <- merge(zoopMets_test, ecoreg, by='UID') |>
              reshape(idvar = c('UID', 'ECO_BIO'), direction = 'wide',
                      v.names = 'RESULT', timevar = 'PARAMETER')

            names(mmiIn) <- gsub("RESULT\\.", "", names(mmiIn))

            testOut <- calcNLA_ZoopMMI(mmiIn, sampID = 'UID', ecoreg = 'ECO_BIO',
                                       totlnind = 'TOTL_NIND')
            varLong <- names(testOut)[names(testOut) %nin% c('UID','ECO_BIO')]
            testOut.long <- reshape(testOut, idvar = c('UID', 'ECO_BIO'), direction = 'long',
                                    varying = varLong,timevar = 'PARAMETER',
                                    v.names = 'RESULT', times = varLong)
            testOut.long <- subset(testOut.long, !is.na(RESULT))

            varLong.1 <- names(zoopMMI_test)[names(zoopMMI_test) %nin% c('UID','ECO_BIO')]
            zoopMMI_test.long <- reshape(zoopMMI_test, idvar = 'UID', direction = 'long',
                                             varying = varLong.1, times = varLong.1, timevar = 'PARAMETER',
                                             v.names = 'RESULT')
            zoopMMI_test.long <- subset(zoopMMI_test.long, !is.na(RESULT))
            compOut <- merge(zoopMMI_test.long, testOut.long, by=c('UID','PARAMETER'))
            compOut.cond <- subset(compOut, PARAMETER=='ZOOP_MMI_COND')
            expect_equal(compOut.cond$RESULT.x,compOut.cond$RESULT.y)
            compOut.res <- subset(compOut,PARAMETER!='ZOOP_MMI_COND')
            compOut.res$RESULT.x <- as.numeric(compOut.res$RESULT.x)
            compOut.res$RESULT.y <- as.numeric(compOut.res$RESULT.y)
            expect_equal(compOut.res$RESULT.x,compOut.res$RESULT.y,tolerance=0.01)

          })
