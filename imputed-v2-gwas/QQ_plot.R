#! /usr/bin/Rscript

## QQ_plot.R

## Author: Daniel P. Howrigan (daniel.howrigan@gmail.com)
## Last Modified: Sept 2017

## Example usage:
## Rscript QQ_plot.R [summary stat file] [output name]

## summary stat file: summary stat .tsv file exported from Hail linreg3
## output: output filename root


## ==== part 1: plotting functions
## ==== part 2: read in summary stats
## ==== part 3: plot summary stats


## ==== part 1: plotting functions

# Stephen Turner
# http://StephenTurner.us/
# http://GettingGeneticsDone.blogspot.com/
# See license at http://gettinggeneticsdone.blogspot.com/p/copyright.html

# Last updated: Tuesday, April19, 2011
# R code for making manhattan plots and QQ plots from plink output files. 
# manhattan() with GWAS data this can take a lot of memory, recommended for use on 64bit machines only, for now. 
# Altnernatively, use bmanhattan() , i.e., base manhattan. uses base graphics. way faster.


## Make a pretty QQ plot of p-values
qq = function(pvector, ...) {
  if (!is.numeric(pvector)) stop("D'oh! P value vector is not numeric.")
  pvector <- pvector[!is.na(pvector) & pvector<1 & pvector>0]
  o = -log10(sort(pvector,decreasing=F))
	e = -log10( ppoints(length(pvector) ))
	plot(e,o,pch=19,cex=1, xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))), xlim=c(0,max(e)), ylim=c(0,max(o)), ...)
	abline(0,1,col="red")
}


## ==== part 2: read in summary stats

## get data.table loaded into R

if(require("data.table")){
    print("data.table is loaded correctly")
} else {
    print("trying to install data.table")
    install.packages("data.table")
    if(require(data.table)){
        print("data.table installed and loaded")
    } else {
        stop("could not install data.table")
    }
}

args <- commandArgs(TRUE)
stats <- args[1] ## summary stat file
output <- args[2] ## output name

gwas <- fread(paste0('zcat ',stats),h=T,stringsAsFactors=F)



## ==== part 3: plot summary stats


## ------ QQ plot
png(paste0(output,'_QQ.png'),height=600,width=600,res=100)

qq(gwas$pval)

lmbx <- qchisq(median(gwas$pval,na.rm=T),1,lower.tail=FALSE)/qchisq(0.50,1,lower.tail=FALSE)

legend('topleft',legend=bquote("Median"~lambda == .(round(lmbx,2))),bg='white',bty='n',cex=1)

dev.off()


