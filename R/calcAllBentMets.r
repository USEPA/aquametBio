#' @export
#'
#' @title Calculate all benthic metrics
#'
#' @description This is a wrapper function that calculates
#' all benthic metrics as used in the National Aquatic
#' Resource Surveys.
#'
#' @param indf A data frame containing, at minimum, the variables
#' specified in the arguments for sampID, dist, ct, and taxa_id
#' @param inTaxa a data frame containing taxonomic information,
#' including variables for PHYLUM, CLASS, ORDER, FAMILY, SUBFAMILY,
#' TRIBE, and GENUS, as well as autecology traits with names that match those
#' in the arguments ffg, habit, and ptv. In addition, there
#' should be a variable with the name in argument taxa_id that matches
#' with all of those in the indf data frame.
#' @param sampID A character vector containing the names of all
#' variables in indf that specify a unique sample. If not specified,
#' the default is \emph{UID}
#' @param dist A string with the name of the distinctness variable,
#' which is assumed to have only values of 0 or 1. If not specified,
#' the default is \emph{IS_DISTINCT}.
#' @param ct A string with the name of the count variable. If not
#' specified, the default is \emph{TOTAL}.
#' @param taxa_id A string with the name of the taxon ID variable
#' in \emph{indf} that matches that in \emph{inTaxa}. The default
#' value is \emph{TAXA_ID}.
#' @param ffg A string with the name of the functional feeding group
#' variable in inTaxa. The default value is \emph{FFG}. Values used
#' in calculations include CF, CG, PR, SH, Sc, representing
#' collector-filterer, collector-gatherer, predator, shredder, and
#' scraper, respectively. Each taxon may have more than
#' one FFG value.
#' @param habit A string with the name of the habit variable in inTaxa.
#' The default value is \emph{HABIT}. Values for habit that are used in
#' calculations include BU, CB, CN, SP, SW, representing burrower,
#' climber, clinger, sprawler, and swimmer, respectively. Each taxon
#' may have more than one value for HABIT.
#' @param ptv A string with the name of the pollution tolerance value
#' variable in inTaxa. The default is \emph{PTV}.
#' @return A data frame containing the variables in sampID and all of
#' the benthic macroinvertebrate metrics as additional variables.
#' @author Karen Blocksom \email{Blocksom.Karen@epa.gov}
#' @examples
#'   data(bentEx)
#'   head(bentEx)
#'   head(bentTaxa_nrsa)
#'   # Calculate metrics for bentIn, using the taxonomy in the count file as is
#'   bentMetrics <- calcAllBentMets(indf=bentEx, inTaxa=bentTaxa_nrsa,
#'                      sampID=c('UID','SAMPLE_TYPE','SAMPLE_CAT'),
#'                      dist='IS_DISTINCT',
#'                      ct='TOTAL',taxa_id='TAXA_ID',
#'                      ffg='FFG',habit='HABIT',ptv='PTV')
#'   head(bentMetrics)
#' @keywords survey
calcAllBentMets <- function(indf,inTaxa, sampID="UID", dist="IS_DISTINCT",
                        ct="TOTAL",taxa_id='TAXA_ID',ffg='FFG',habit='HABIT',ptv='PTV'){

  # Make sure all taxa match to taxalist and send error if not
  checkTaxa <- indf[indf$TAXA_ID %in% setdiff(indf$TAXA_ID, inTaxa$TAXA_ID),]
  #checkTaxa <- dplyr::anti_join(indf,inTaxa,by='TAXA_ID')
  if(nrow(checkTaxa)>0){
    return(print('Taxa in counts that do not have matches in taxalist! Cannot continue.'))
  }

  if('NON_TARGET' %in% names(inTaxa)){
    inTaxa <- subset(inTaxa, is.na(NON_TARGET) | NON_TARGET == "" |NON_TARGET=='N')
  }

  ctVars <- c(sampID,dist,ct,taxa_id)
  if(any(ctVars %nin% names(indf))){
    msgTraits <- which(ctVars %nin% names(indf))
    print(paste("Missing variables in input count data frame:",paste(names(indf)[msgTraits],collapse=',')))
    return(NULL)
  }

  # Taxonomy and traits checks
  necTraits <- c('PHYLUM','CLASS','ORDER','FAMILY','TRIBE','SUBFAMILY','GENUS'
                 ,ffg,habit,ptv)
  if(any(necTraits %nin% names(inTaxa))){
    msgTraits <- which(necTraits %nin% names(inTaxa))
    return(paste("Some of the traits are missing from the taxa list. The following are \nrequired for metric calculations to run:\n", necTraits[msgTraits], "\n"))
  }

  inTaxa <- subset(inTaxa,select=names(inTaxa) %in% c('TAXA_ID','PHYLUM','CLASS','ORDER','FAMILY'
                                                      ,'TRIBE','SUBFAMILY','GENUS',ffg,habit,ptv))

  indf[,c(ct,taxa_id,dist)] <- lapply(indf[,c(ct,taxa_id,dist)],as.numeric)
  inTaxa[,c(ptv,taxa_id)] <- lapply(inTaxa[,c(ptv,taxa_id)],as.numeric)


  taxMet <- calcBentTaxMets(indf,inTaxa,sampID,dist,ct,taxa_id)
  tax.1 <- reshape(taxMet, idvar = sampID, direction = 'long', varying = c('AMPHNTAX','AMPHPIND'
                                                                           ,'AMPHPTAX','CHIRNTAX','CHIRPIND'
                                                                           ,'CHIRPTAX','CRUSNTAX','CRUSPIND','CRUSPTAX'
                                                                           ,'DIPTNTAX','DIPTPIND','DIPTPTAX','EPHENTAX'
                                                                           ,'EPHEPIND','EPHEPTAX','EPOTNTAX','EPOTPIND'
                                                                           ,'EPOTPTAX','EPT_NTAX','EPT_PIND'
                                                                           ,'EPT_PTAX','HEMINTAX','HEMIPIND','HEMIPTAX'
                                                                           ,'MITENTAX','MITEPIND','MITEPTAX'
                                                                           ,'MOLLNTAX','MOLLPIND','MOLLPTAX','NOINNTAX','NOINPIND','NOINPTAX'
                                                                           ,'ODONNTAX','ODONPIND','ODONPTAX','OLLENTAX','OLLEPIND'
                                                                           ,'OLLEPTAX','ORTHNTAX','ORTHPIND','ORTHPTAX'
                                                                           ,'PLECNTAX','PLECPIND','PLECPTAX'
                                                                           ,'TANYNTAX','TANYPIND'
                                                                           ,'TANYPTAX','TRICNTAX','TRICPIND','TRICPTAX'
                                                                           ,'TUBINAIDNTAX','TUBINAIDPIND','TUBINAIDPTAX','ORTHCHIRPIND'),
                   v.names = 'value', timevar = 'variable', times = c('AMPHNTAX','AMPHPIND'
                                                                      ,'AMPHPTAX','CHIRNTAX','CHIRPIND'
                                                                      ,'CHIRPTAX','CRUSNTAX','CRUSPIND','CRUSPTAX'
                                                                      ,'DIPTNTAX','DIPTPIND','DIPTPTAX','EPHENTAX'
                                                                      ,'EPHEPIND','EPHEPTAX','EPOTNTAX','EPOTPIND'
                                                                      ,'EPOTPTAX','EPT_NTAX','EPT_PIND'
                                                                      ,'EPT_PTAX','HEMINTAX','HEMIPIND','HEMIPTAX'
                                                                      ,'MITENTAX','MITEPIND','MITEPTAX'
                                                                      ,'MOLLNTAX','MOLLPIND','MOLLPTAX','NOINNTAX','NOINPIND','NOINPTAX'
                                                                      ,'ODONNTAX','ODONPIND','ODONPTAX','OLLENTAX','OLLEPIND'
                                                                      ,'OLLEPTAX','ORTHNTAX','ORTHPIND','ORTHPTAX'
                                                                      ,'PLECNTAX','PLECPIND','PLECPTAX'
                                                                      ,'TANYNTAX','TANYPIND'
                                                                      ,'TANYPTAX','TRICNTAX','TRICPIND','TRICPTAX'
                                                                      ,'TUBINAIDNTAX','TUBINAIDPIND','TUBINAIDPTAX','ORTHCHIRPIND'))
  # tax.1 <- data.table::melt(taxMet,id.vars=sampID)

  ffgMet <- calcBentFFGmets(indf,inTaxa,sampID,dist,ct,taxa_id,ffg)
  ffg.1 <- reshape(ffgMet, idvar = sampID, direction = 'long', varying = c('COFINTAX','COFIPIND','COFIPTAX'
                                                                           ,'COFITRICNTAX','COFITRICPIND','COFITRICPTAX'
                                                                           ,'COGANTAX','COGAPIND','COGAPTAX'
                                                                           ,'PREDNTAX','PREDPIND','PREDPTAX'
                                                                           ,'SCRPNTAX','SCRPPIND','SCRPPTAX'
                                                                           ,'SHRDNTAX','SHRDPIND','SHRDPTAX'),
                   v.names = 'value', timevar = 'variable', times = c('COFINTAX','COFIPIND','COFIPTAX'
                                                                      ,'COFITRICNTAX','COFITRICPIND','COFITRICPTAX'
                                                                      ,'COGANTAX','COGAPIND','COGAPTAX'
                                                                      ,'PREDNTAX','PREDPIND','PREDPTAX'
                                                                      ,'SCRPNTAX','SCRPPIND','SCRPPTAX'
                                                                      ,'SHRDNTAX','SHRDPIND','SHRDPTAX'))
  # ffg.1 <- data.table::melt(ffgMet,id.vars=sampID)

  habitMet <- calcBentHabitMets(indf,inTaxa,sampID,dist,ct,taxa_id,habit)
  habit.1 <- reshape(habitMet, idvar = sampID, direction = 'long', varying = c('BURRNTAX','BURRPIND','BURRPTAX'
                                                                               ,'CLMBNTAX','CLMBPIND','CLMBPTAX'
                                                                               ,'CLNGNTAX','CLNGPIND','CLNGPTAX'
                                                                               ,'SPWLNTAX','SPWLPIND','SPWLPTAX'
                                                                               ,'SWIMNTAX','SWIMPIND','SWIMPTAX'),
                     v.names = 'value', timevar = 'variable', times = c('BURRNTAX','BURRPIND','BURRPTAX'
                                                                        ,'CLMBNTAX','CLMBPIND','CLMBPTAX'
                                                                        ,'CLNGNTAX','CLNGPIND','CLNGPTAX'
                                                                        ,'SPWLNTAX','SPWLPIND','SPWLPTAX'
                                                                        ,'SWIMNTAX','SWIMPIND','SWIMPTAX'))
  # habit.1 <- data.table::melt(habitMet,id.vars=sampID)

  tolMet <- calcBentTolMets(indf,inTaxa,sampID,dist,ct,taxa_id,ptv)
  tol.1 <- reshape(tolMet, idvar = sampID, direction = 'long', varying = c('FACLNTAX','FACLPIND','FACLPTAX'
                                                                           ,'INTLNTAX','INTLPIND','INTLPTAX'
                                                                           ,'NTOLNTAX','NTOLPIND','NTOLPTAX'
                                                                           ,'STOLNTAX','STOLPIND','STOLPTAX'
                                                                           ,'TL01NTAX','TL01PIND','TL01PTAX'
                                                                           ,'TL23NTAX','TL23PIND','TL23PTAX'
                                                                           ,'TL45NTAX','TL45PIND','TL45PTAX'
                                                                           ,'TL67NTAX','TL67PIND','TL67PTAX'
                                                                           ,'TOLRNTAX','TOLRPIND','TOLRPTAX'
                                                                           ,'WTD_TV'),
                   v.names = 'value', timevar = 'variable', times = c('FACLNTAX','FACLPIND','FACLPTAX'
                                                                      ,'INTLNTAX','INTLPIND','INTLPTAX'
                                                                      ,'NTOLNTAX','NTOLPIND','NTOLPTAX'
                                                                      ,'STOLNTAX','STOLPIND','STOLPTAX'
                                                                      ,'TL01NTAX','TL01PIND','TL01PTAX'
                                                                      ,'TL23NTAX','TL23PIND','TL23PTAX'
                                                                      ,'TL45NTAX','TL45PIND','TL45PTAX'
                                                                      ,'TL67NTAX','TL67PIND','TL67PTAX'
                                                                      ,'TOLRNTAX','TOLRPIND','TOLRPTAX'
                                                                      ,'WTD_TV'))
  # tol.1 <- data.table::melt(tolMet,id.vars=sampID)

  domMet <- calcBentDominMets(indf,inTaxa,sampID,dist,ct,taxa_id)
  # Still working on this below - do I know the set variable names ahead of time?
  dom.1 <- reshape(domMet, idvar = sampID, direction='long', varying=c('HPRIME','DOM1PIND','DOM3PIND','DOM5PIND',
                                                                       'CHIRDOM1PIND','CHIRDOM3PIND','CHIRDOM5PIND'),
                   v.names='value', timevar='variable', times=c('HPRIME','DOM1PIND','DOM3PIND','DOM5PIND','CHIRDOM1PIND',
                                                                'CHIRDOM3PIND','CHIRDOM5PIND'))
  # dom.1 <- data.table::melt(domMet,id.vars=sampID)

  names(indf)[names(indf)==ct] <- 'FINAL_CT'
  names(indf)[names(indf)==dist] <- 'IS_DISTINCT'

  rhs <- paste(sampID, collapse='+')
  form <- paste('cbind(FINAL_CT, IS_DISTINCT)', rhs, sep='~')
  totals <- aggregate(formula(form), data=indf, FUN = function(x) sum=sum(x))

  names(totals)[names(totals)=='FINAL_CT'] <- 'TOTLNIND'
  names(totals)[names(totals)=='IS_DISTINCT'] <- 'TOTLNTAX'

  # Now melt the above using variable and value?
  totals.long <- reshape(totals, idvar = sampID, direction='long', varying=c('TOTLNIND','TOTLNTAX'),v.names='value',timevar='variable',times=c('TOTLNIND','TOTLNTAX'))

  # totals <- plyr::ddply(indf, sampID, summarise, TOTLNIND=sum(FINAL_CT),
  #                                  TOTLNTAX=sum(IS_DISTINCT)) %>%
  #   melt(id.vars=sampID)

  mets <- rbind(tax.1, ffg.1, habit.1, tol.1, dom.1,totals.long)

  # Finally, we can recast the metrics df into wide format for output
  lside <- paste(paste(sampID,collapse='+'),sep='+')
  form <- paste(lside,'~variable',sep='')
  metOut <- reshape(mets, direction = 'wide', idvar = sampID,  timevar = 'variable')
  names(metOut) <- gsub("value\\.", "", names(metOut))
  # metOut <- data.table::dcast(mets,eval(form),value.var='value')

  return(metOut)
}
