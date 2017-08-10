#!/bin/sh

BGEN_BUCKET="gs://fc-9a7c5487-04c9-4182-b3ec-13de7f6b409b/"
WORKING_BUCKET="gs://ukbb_association/"

gsutil -m cp ${BGEN_BUCKET}/imputed/ukb_mfi_chr*_v2.txt ./

for i in {1..22}; do
	awk -v chr=$i 'BEGIN {FS="\t"; OFS="\t"} {print chr,$0}' "ukb_mfi_chr${i}_v2.txt" >> ukb_mfi_v2.tsv
done

gzip ukb_mfi_v2.tsv

gsutil cp ukb_mfi_v2.tsv.gz ${WORKING_BUCKET}
