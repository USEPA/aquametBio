###############################################################################
# metCalcFxns.r
#
# Contains functions used in metric calculation programs for fish and bugs
#
# Programmers: Karen Blocksom
#              Tom Kincaid
# Date: June 5, 2013
# Function Revisions:
#   06/05/13 kab: Created based on MetricCalcFxns.R.
#   08/09/13 kab: Updated to use only distinct taxa in Dominance() calculations,
#            following A. Herlihy calculations.
#   03/04/14 kab: Updated to keep old version of tolerance index for fish as
#            tolindexFish()
#   03/05/14 tmk: Removed call to the require() function.
###############################################################################


# This function calculates % dominant organisms metric
# df <-- input data frame, containing UID as variable identifying unique samples
# topN <-- number specifying the top number of species to include in calculation

Dominance<-function(df, topN=1){
 rr <- subset(df,IS_DISTINCT==1)
 totals <- aggregate(list(TOTSUM = rr$TOTAL), by = rr[c('SAMPID')], FUN = sum )

 ss <- merge(rr, totals, by = 'SAMPID')

	tt <- stats::aggregate(list(domN=ss$TOTAL)
			,list(SAMPID=ss$SAMPID)
			,function(x){
				sum(x[order(x,decreasing=TRUE)[1:topN]]
				)
			}
	)
	uu <- merge(tt,unique(ss[,c('SAMPID','TOTSUM')]),by="SAMPID")
	uu$dompind <- with(uu, round(domN/TOTSUM*100,2))
	uu <- subset(uu,select=c('SAMPID','dompind'))
	names(uu)[names(uu) == 'dompind'] <- paste('DOM', topN, 'PIND', sep = '')

	return(uu)
}


# Function calculates Shannon Diversity metric
# indata <-- data frame containing unique sample identifier (UID), IS_DISTINCT numeric variable indicating taxonomic distinctness as 0 or 1,
# 		and TOTAL as the count variable

ShanDiversity <- function(indata){
	rr <- subset(indata,IS_DISTINCT==1)
	ss <- aggregate(x=list(TOTSUM = rr$TOTAL), by = rr[c('SAMPID')], FUN = sum)

	rr <- merge(rr, ss, by = 'SAMPID')
	rr$CALC <- with(rr,(TOTAL/TOTSUM)*(log(TOTAL/TOTSUM)))

	tt <- aggregate(x = list(HPRIME = rr$CALC), by = rr[c('SAMPID')], FUN = function(x){ round(-1*sum(x), 2)})

	return(tt)
}


# Function to calculate weighted tolerance value index (e.g., HBI) as the sum of the proportion of a taxon multiplied by the tolerance value
#		for that taxon. All taxa included in proportion calculations, not just those with PTVs.
# indata <-- input data frame with TAXA_ID, TOTAL, and UID variables
# taxalist <-- taxalist with TAXA_ID and tolerance values (TVs) as RESULT and PARAMETER as either PTV or TOL_VAL

tolindex<-function(indata, taxalist){
	tv_taxa <- taxalist[taxalist$TRAIT %in% c('PTV','TOL_VAL'),]
	rr <- aggregate(x = list(SUMCT = indata$TOTAL), by = indata[c('SAMPID')], FUN = sum)
	#This allows us to sum across only those taxa with TVs
	tv_cts <- merge(indata,tv_taxa,by="TAXA_ID")
	# Redo this to match order of operations of original code - get total separately
	tv_cts$CALC <- with(tv_cts, TOTAL*as.numeric(value))
	tvOut <- aggregate(x = list(SUMCALC = tv_cts$CALC), by = tv_cts[c('SAMPID')], FUN = sum)
	outTV <- merge(tvOut, rr, by = c('SAMPID'))
	outTV$WTD_TV <- with(outTV, round(SUMCALC/SUMCT, 2))
	return(outTV)
}


# Original version: Only taxa with TVs are included in the proportion
# calculations.
tolindexFish<-function(indata, taxalist){
	tv_taxa <- taxalist[,c('TAXA_ID','TOL_VAL')]
	# This allows us to sum across only those taxa with TVs
	tv_cts <- merge(indata,tv_taxa,by="TAXA_ID")

	rr <- subset(tv_cts, !is.na(TOTAL) & !is.na(TOL_VAL))
	totsum <- aggregate(x = list(TOTSUM = rr$TOTAL), by = rr[c('SAMPID')], FUN = function(x){ sum(x, na.rm=TRUE) })

	ss <- merge(rr, totsum, by = 'SAMPID')
	ss$CALC <- with(ss, TOTAL*as.numeric(TOL_VAL)/TOTSUM)

	outTV <- aggregate(x = list(WTD_TV = ss$CALC), by = ss[c('SAMPID')], FUN = function(x){ round(sum(x), 2)})
	return(outTV)
}
