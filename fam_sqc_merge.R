#!/usr/bin/Rscript

## fam_sqc_merge.R

## Arguments to read in:

# Rscript fam_sqc_merge.R \
# [UKBB qc file] \
# [UKBB header file] \
# [application specific .fam file] \
# [output file name]

## Example usage:

# Rscript fam_sqc_merge.R \
# ukb_sqc_v2.txt \
# ukb_sqc_v2_header.txt \
# ukb1859_cal_chr22_v2_s488374.fam \
# ukb1859_qc.txt

args <- commandArgs(TRUE)
sqc_file <- args[1]
hdr_file <- args[2]
fam_file <- args[3]
out_file <- args[4]

library(data.table)
sqc <- fread(sqc_file,stringsAsFactors=F)
sqc <- data.frame(sqc)
sqc <- sqc[,3:ncol(sqc)]
fam <- fread(fam_file)
fam <- data.frame(fam)
eid <- fam$V1

sqc2 <- cbind.data.frame(eid,sqc)

hd <- read.table(hdr_file)
## add sample names
hd <- c(c('eid'),as.character(hd$V1)[3:nrow(hd)])

if (ncol(sqc2)==length(hd)) {
	names(sqc2) <- hd
	write.table(sqc2,out_file,col=T,row=F,quo=F,sep='\t')
} 

if (ncol(sqc2)!=length(hd)) {
	cat('header and qc file are not concordant\n')
} 

