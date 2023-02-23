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
