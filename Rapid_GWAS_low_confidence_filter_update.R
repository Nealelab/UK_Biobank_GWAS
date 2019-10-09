## R script: Rapid_GWAS_low_confidence_filter_update.R

#! bin/Rscript

## author: Daniel P. Howrigan
## email: howrigan@broadinstitute.org


## USAGE:
# Rscript Rapid_GWAS_low_confidence_filter_update.R [summary file] [summary statistic directory] [output directory]


## EXAMPLE USAGE:
# Rscript Rapid_GWAS_low_confidence_filter_update.R \
# GWAS_list_low_confidence_filter_update.txt.gz \
# rapid_gwas_sumstats/ \
# updated_rapid_gwas_sumstats/



## get data.table loaded into R
## NOTE: currently using "data.table 1.12.2" 

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



## read in arguments
args <- commandArgs(TRUE)
phe_file_name <- args[1] ## GWAS filename list
in_dir <- args[2] ## input directory
out_dir <- args[2] ## output directory

# confirm/create output dir
system(paste0('mkdir -p ',out_dir))


## read in file

phe_file <- fread(cmd=paste0('zcat ',phe_file_name),h=T,quote='',stringsAsFactors=F) 

## create variables to collect
gwas_file <- NA
prop_low_confidence_in <- NA
mean_AC_low_confidence_in <- NA
smallest_cat_size_used <- NA
smallest_cat_size_PHESANT <- NA
prop_low_confidence_out <- NA
mean_AC_low_confidence_out <- NA


## Loop through phenotypes

for (i in 1:nrow(phe_file)) {

## read in GWAS sumstat
gwas <- fread(cmd=paste0('zcat ',in_dir,phe_file$additive_tsvs_list_name[i]),h=T,stringsAsFactors=F)

gwas_file[i] <- phe_file$additive_tsvs_list_name[i]
# prop_low_confidence_in[i] <- sum(gwas$low_confidence_variant=='true') / nrow(gwas)
# mean_AC_low_confidence_in[i] <- mean(gwas$AC[gwas$low_confidence_variant=='true'])
prop_low_confidence_in[i] <- sum(gwas$low_confidence_variant==T) / nrow(gwas)
mean_AC_low_confidence_in[i] <- mean(gwas$AC[gwas$low_confidence_variant==T])

## calculate smallest_cat_size_used with 1:693731:A:G
snp <- subset(gwas,gwas$variant=='1:693731:A:G')
smallest_cat_size_used[i] <- round(snp$expected_min_category_minor_AC / ((snp$minor_AF) * 2),0)

smallest_cat_size_PHESANT[i] <- phe_file$min_category[i]

## create new expected AC and FILTER
expected_min_category_minor_AC <- gwas$minor_AF * phe_file$min_category[i] * 2
# low_confidence_variant <- rep('false',nrow(gwas))
# low_confidence_variant[expected_min_category_minor_AC < 25] <- 'true'
# low_confidence_variant[gwas$minor_AF < 0.001] <- 'true'

# prop_low_confidence_out[i] <- sum(low_confidence_variant=='true') / nrow(gwas)
# mean_AC_low_confidence_out[i] <- mean(gwas$AC[low_confidence_variant=='true'])

low_confidence_variant <- rep(FALSE,nrow(gwas))
low_confidence_variant[expected_min_category_minor_AC < 25] <- TRUE
low_confidence_variant[gwas$minor_AF < 0.001] <- TRUE

prop_low_confidence_out[i] <- sum(low_confidence_variant==TRUE) / nrow(gwas)
mean_AC_low_confidence_out[i] <- mean(gwas$AC[low_confidence_variant==TRUE])


## overwrite AC and FILTER
gwas$expected_min_category_minor_AC <- expected_min_category_minor_AC
gwas$low_confidence_variant <- low_confidence_variant


## write updated GWAS file
nm <- phe_file$additive_tsvs_list_name[i]
nm2 <- substr(phe_file$additive_tsvs_list_name[i],start=1,stop=nchar(phe_file$additive_tsvs_list_name[i])-4)
write.table(gwas,paste0(out_dir,nm2),col=T,row=F,quo=F,sep='\t')

## gzip file
system(paste0("gzip ",out_dir,nm2))

## print to stdout
cat(paste0(i,' of ',nrow(phe_file),' complete'),'\n')
cat(paste0('gwas_file = ',gwas_file[i]),'\n')
cat(paste0('prop. low conf. variants = ',round(prop_low_confidence_in[i],2)),'\n')
cat(paste0('mean AC low conf. variants = ',round(mean_AC_low_confidence_in[i],1)),'\n')
cat(paste0('smallest cat used = ',smallest_cat_size_used[i]),'\n')
cat(paste0('smallest cat in PHESANT = ',smallest_cat_size_PHESANT[i]),'\n')
cat(paste0('prop. low conf. variants after update = ',round(prop_low_confidence_out[i],2)),'\n')
cat(paste0('mean AC low conf. variants after update = ',round(mean_AC_low_confidence_out[i],1)),'\n')


} ## END of i LooP


summary_file <- cbind.data.frame(
gwas_file,
prop_low_confidence_in,
mean_AC_low_confidence_in,
smallest_cat_size_used,
smallest_cat_size_PHESANT,
prop_low_confidence_out,
mean_AC_low_confidence_out)

## write to file
nm <- substr(phe_file_name,start=1,stop=nchar(phe_file_name)-4)
write.table(summary_file,paste0(nm,'_summary.tsv'),col=T,row=F,quo=F,sep='\t')

## ===== END of R script


