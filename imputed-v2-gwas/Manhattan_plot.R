#! /usr/bin/Rscript

## Manhattan_plot.R

## Author: Daniel P. Howrigan (daniel.howrigan@gmail.com)
## Last Modified: Sept 2017

## Example usage:
## Rscript Manhattan_plot.R [summary stat file] [output name]

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


## Code edited by Daniel P. Howrigan


# manhattan plot using base graphics
manhattan <- function(dataframe, colors=c("gray10", "gray50"), ymax="max", limitchromosomes=1:23, suggestiveline=-log10(1e-5), genomewideline=-log10(5e-8), annotate=NULL, ...) {

    d=dataframe
    if (!("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, and P")
    
    if (any(limitchromosomes)) d=d[d$CHR %in% limitchromosomes, ]
    d=subset(na.omit(d[order(d$CHR, d$BP), ]), (P>0 & P<=1)) # remove na's, sort, and keep only 0<P<=1
    d$logp = -log10(d$P)
    d$pos=NA
    ticks=NULL
    lastbase=0
    colors <- rep(colors,max(d$CHR))[1:max(d$CHR)]
    if (ymax=="max") ymax<-ceiling(max(d$logp))
    if (ymax<8) ymax<-8
    
    numchroms=length(unique(d$CHR))
    if (numchroms==1) {
        d$pos=d$BP
        ticks=floor(length(d$pos))/2+1
    } else {
        for (i in unique(d$CHR)) {
          if (i==1) {
    			d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP
    		} else {
    			lastbase=lastbase+tail(subset(d,CHR==i-1)$BP, 1)
    			d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP+lastbase
    		}
    		ticks=c(ticks, d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
    	}
    }
    
    if (numchroms==1) {
        with(d, plot(pos, logp, ylim=c(0,ymax), ylab=expression(-log[10](italic(p))), xlab=paste("Chromosome",unique(d$CHR),"position"), ...))
    }	else {
        with(d, plot(pos, logp, ylim=c(0,ymax), ylab=expression(-log[10](italic(p))), xlab="Chromosome", xaxt="n", type="n", ...))
        axis(1, at=ticks, lab=unique(d$CHR), ...)
        icol=1
        for (i in unique(d$CHR)) {
            with(d[d$CHR==i, ],points(pos, logp, col=colors[icol], ...))
            icol=icol+1
    	}
    }
    
    if (!is.null(annotate)) {
        d.annotate=d[which(d$SNP %in% annotate), ]
        with(d.annotate, points(pos, logp, col="green3", ...)) 
    }
    
    if (suggestiveline) abline(h=suggestiveline, col="blue")
    if (genomewideline) abline(h=genomewideline, col="red")
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

## subset to significant hits
## and desired columns
awkcmd = 'awk \'NR==1 { for (i=1; i<=NF; i++) {f[$i] = i} } { if(NR==1 || $f["pval"] < 0.001) {print $f["variant"], $f["pval"]} }\''
gwas <- fread(paste("gunzip -c",stats,"|",awkcmd),h=T,stringsAsFactors=F)


## ==== part 3: plot summary stats


## ------ Manhattan plot

## Just keep CHR, BP, and P
tmp <- strsplit(gwas$variant,':',fixed=T)
CHR <- as.numeric(sapply(tmp,'[',1))
BP <- as.numeric(sapply(tmp,'[',2))
P <- gwas$pval
P[is.na(P)] <- 1

gwas3 <- cbind.data.frame(CHR,BP,P)

png(paste0(output,'_Manhattan.png'),height=600,width=1200,res=100,pointsize=3)

manhattan(gwas3,pch=16,cex=3,colors=c("lightskyblue2","midnightblue"),suggestiveline=FALSE,cex.axis=3,main='')

dev.off()


