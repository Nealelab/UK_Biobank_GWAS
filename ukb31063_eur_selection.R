#! /usr/bin/env Rscript

# setup 
setwd("/home/rstudio/rkw/data/ukb31063/eur_selection/")
require(data.table,lib.loc="/home/rstudio/R/x86_64-pc-linux-gnu-library/3.4")
require(RColorBrewer,lib.loc="/home/rstudio/R/x86_64-pc-linux-gnu-library/3.4")
ells <- 6
sds_white <- 4.5
sds_brit <- 7


#####
# load data
#####

# load QC data
qc <- fread("ukb31063_sample_qc.tsv", sep='\t', header=T, stringsAsFactors=F, data.table=F)

# filter QC criteria
qc_pass <- (qc$het.missing.outliers==0 &
              qc$excluded.from.kinship.inference==0 &
              qc$excess.relatives==0 &
              qc$used.in.pca.calculation==1)

qc2 <- qc[qc_pass,]


# select phenotypes to keep
ph_cols <- c(
  userId="userId",
  country="x1647_0_0", # country of birth in UK
  country2="x20115_0_0", # country of birth outside UK
  ethnic="x21000_0_0" # self-report ethnicity
) 

# read selected phenotypes
phens <- fread("neale_lab_parsed.tsv", 
               sep='\t', 
               header=T, 
               stringsAsFactors=F, 
               data.table=F,
               select=unname(ph_cols))
names(phens) <- names(ph_cols)

#merge info
df <- merge(phens, qc2, by.y="iid", by.x="userId")


# get self-report whites
df$white <- (df$ethnic %in% c(1,1001,1002,1003))
df$ethnic_miss <- (df$ethnic %in% c(-3,-1,NA))


####
# get white european selection
####

# get mean and SD of each PC among 
# the curated white British sample
# and self-reported whites
pc_nams <- paste("PC",1:40,sep="")
mm_white <- colMeans(df[df$white==1,pc_nams])
ss_white <- apply(df[df$white==1,pc_nams],2,sd)
mm_brit <- colMeans(df[df$in.white.British.ancestry.subset==1,pc_nams])
ss_brit <- apply(df[df$in.white.British.ancestry.subset==1,pc_nams],2,sd)

# draw ellipses
dd_white <- rep(0,nrow(df))
dd_brit <- rep(0,nrow(df))
for(i in 1:ells){
  dd_white <- dd_white + (df[,pc_nams[i]]-mm_white[i])^2/(ss_white[i]^2)
  dd_brit <- dd_brit + (df[,pc_nams[i]]-mm_brit[i])^2/(ss_brit[i]^2)
}

# make selection
# intersection of:
# - curated white british ellipse
# - self-reported white ellipse
# - self-reported white
# df$eur_select <- (dd_white < sds_white^2) & (dd_brit < sds_brit^2) & (df$white | df$ethnic_miss)
df$eur_select <- (dd_brit < sds_brit^2) & (df$white | df$ethnic_miss)

# save selection
write.table(df[,c("userId","eur_select")],file="ukb31063_eur_samples.tsv",sep='\t',col.names=T,row.names=F)


####
# plot
####

# colors
set <- rep(1,nrow(df))
set[df$eur_select] <- 3
set[df$in.white.British.ancestry.subset==1] <- 2
cols <- c("gray80",brewer.pal("Dark2",n=3)[1],brewer.pal("Dark2",n=3)[2])

# samp=sample(nrow(df),10000,replace=F)
samp <- 1:nrow(df)
png("ukb31063_eur_selection_pca.png",width=18,height=18,res=300,units="in")
  pairs(df[samp,c("PC1","PC2","PC3","PC4","PC5","PC6")],col=cols[set[samp]],cex=.1)
dev.off()

###
# summaries
###

print("versus white British")
table(df$in.white.British.ancestry.subset,df$eur_select,useNA="al")

print("versus reported ethnicity")
table(df$ethnic,set,useNA="al")

print("versus reported country of birth in UK")
table(df$country,set,useNA="if")

print("versus reported country of birth outside UK")
table(df$country2,set,useNA="if")

# eof