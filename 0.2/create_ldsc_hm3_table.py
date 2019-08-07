#! /usr/bin/env python

#############

# Gets set of SNPs from UKB GWAS to export for ldsc
# - UKB MAF > 1%
# - UKB INFO > 0.9
# - autosomal SNP
# - biallelic (in HM3, in UKB MAF/INFO filtered set)
# - in HapMap3 (match location + ref + alt + rsid)
# - not in MHC

# this version using Hail 0.2 
# and variant QC from round2 GWAS

#############



####
# setup
####

# re-run everything?
force_all = True

# count variants at each step?
verbose = True

### input

# source of HapMap3 sites file (path only)
HM3_FTP = 'ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/'

# name of original HapMap3 sites file
HM3_NAME = 'hapmap_3.3_b37_pop_stratified_af.vcf.gz'

# UKBB variant QC in the full sample (for INFO scores)
GWAS_mfi = 'gs://ukb31063/ukb31063.neale_gwas_variants.imputed_v3.mfi.ht'

# UKBB variant QC in the GWAS samples (for MAF)
# Note: using saved QC metrics that match values before last round of sample withdrawals
GWAS_qc = 'gs://ukb-mega-gwas-results/round2/annotations/variants.tsv.bgz'

# SNPs passing QC to be included in the UKBB GWAS
GWAS_snps = 'gs://ukb31063/ukb31063.neale_gwas_variants.ht'


### output

# where to save files
BUCKET = 'gs://ukbb-ldsc-dev/ukb_hm3_snplist/'

# name to save HapMap3 sites as matrixtable
HM3_vds = 'hm3.r3.b37.mt'

# name to save autosomal, biallelic SNPs (also has parsed population allele freqs)
HM3_bi_af = 'hm3.r3.b37.auto_bi_af.ht'

# name to save HapMap3 SNPs passing ldsc QC in UKBB GWAS sample
# - will save both Hail table (.ht) and flat file (.tsv)
HM3_qcpos_stem = 'hm3.r3.b37.auto_bi_af.ukbb_gwas_qcpos.no_mhc'




print("Preparing packages, etc...")

# load packages
import hail as hl
import subprocess
hl.init()


# utility to check if file exists on cloud
def gs_missing(filepath):
    
    # strip tailing '/'
    if str(filepath).endswith('/'):
        tmp = str(filepath)[:-1]
        filepath = tmp
    
    # check for file
    stat = subprocess.call(['gsutil','-q','stat',str(filepath)])
    
    # check for folder
    if stat:
        stat = subprocess.call(['gsutil','-q','stat',str(filepath)+'/'])
    
    return stat


###
# Load hapmap3 sites vcf to hail table
#   hm3 file is from Broad resource bundle:
#   wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/hapmap_3.3_b37_pop_stratified_af.vcf.gz
###

if gs_missing(BUCKET + HM3_NAME) or force_all:
    
    print("Downloading HapMap3 vcf...")
    
    # download from Broad
    subprocess.call(['wget',
                    '--quiet', 
                    '-P', 
                    '/home/hapmap/',
                    HM3_FTP + HM3_NAME])
    
    # upload to cloud
    subprocess.call(['gsutil', 
                     'cp', 
                     '/home/hapmap/' + HM3_NAME,
                     BUCKET + HM3_NAME])
    

if gs_missing(BUCKET + HM3_vds) or force_all: 
    
    print("Creating Hail table of HapMap3 sites...")
        
    (hl
        .import_vcf(BUCKET + HM3_NAME, 
                    force_bgz=True, 
                    reference_genome='GRCh37', 
                    find_replace=("##INFO=<ID=AC,Number=1,","##INFO=<ID=AC,Number=A,"),
                    skip_invalid_loci=False, 
                    min_partitions=500)
        .write(BUCKET + HM3_vds, overwrite=True)
    )



###
# Create filtered keytable with autosomal, biallelic HM3 snps
# with MAF > .01, INFO > 0.9 and passing QC in UKBB GWAS analysis set
###
if gs_missing(BUCKET + HM3_bi_af) or force_all:
    
    ###
    # Create filtered keytable with only autosomal, biallelic snps
    # - also parse to keep allele freqs by population
    print("Creating keytable of autosomal, biallelic HM3 SNPs with population allele freqs...")
    ###
    
    # read full HM3 sites
    hm3_site = hl.read_matrix_table(BUCKET + HM3_vds).rows()
    
    # filter to autosomal, biallelic snps and extract allelic freqs
    hm3_keep = hm3_site.filter(hm3_site.locus.in_autosome() & 
                               (hl.len(hm3_site.alleles) == 2) & 
                               hl.is_snp(hm3_site.alleles[0],hm3_site.alleles[1]))
    
    hm3_keep = hm3_keep.annotate(ASW_AF = 1.0 - hl.float(hm3_keep.info['ASW'][0].split("=")[1]), 
                                 CEU_AF = 1.0 - hl.float(hm3_keep.info['CEU'][0].split("=")[1]),
                                 CHB_AF = 1.0 - hl.float(hm3_keep.info['CHB'][0].split("=")[1]),
                                 CHD_AF = 1.0 - hl.float(hm3_keep.info['CHD'][0].split("=")[1]),
                                 CHS_AF = 1.0 - hl.float(hm3_keep.info['CHS'][0].split("=")[1]),
                                 CLM_AF = 1.0 - hl.float(hm3_keep.info['CLM'][0].split("=")[1]),
                                 FIN_AF = 1.0 - hl.float(hm3_keep.info['FIN'][0].split("=")[1]),
                                 GBR_AF = 1.0 - hl.float(hm3_keep.info['GBR'][0].split("=")[1]),
                                 GIH_AF = 1.0 - hl.float(hm3_keep.info['GIH'][0].split("=")[1]),
                                 IBS_AF = 1.0 - hl.float(hm3_keep.info['IBS'][0].split("=")[1]),
                                 JPT_AF = 1.0 - hl.float(hm3_keep.info['JPT'][0].split("=")[1]),
                                 LWK_AF = 1.0 - hl.float(hm3_keep.info['LWK'][0].split("=")[1]),
                                 MKK_AF = 1.0 - hl.float(hm3_keep.info['MKK'][0].split("=")[1]),
                                 MXL_AF = 1.0 - hl.float(hm3_keep.info['MXL'][0].split("=")[1]),
                                 PUR_AF = 1.0 - hl.float(hm3_keep.info['PUR'][0].split("=")[1]),
                                 TSI_AF = 1.0 - hl.float(hm3_keep.info['TSI'][0].split("=")[1]),
                                 YRI_AF = 1.0 - hl.float(hm3_keep.info['YRI'][0].split("=")[1]))


    # save keytable
    hm3_keep.write(BUCKET + HM3_bi_af, overwrite=True)


###
# Create filtered keytable with autosomal, biallelic HM3 snps
# with MAF > .01, INFO > 0.9 and passing QC in UKBB GWAS analysis set
#
# Also provides both UKBB and HM3 formatting of variant
###

if gs_missing(BUCKET + HM3_qcpos_stem + '.ht') or force_all:
    
    print("Creating Hail table of HM3 SNPs passing UKBB GWAS QC...")


    # get list of SNPs to be used in GWAS
    # filter here: autosomes, no indels, no MHC (filter early for efficiency)
    ukb_snps = hl.read_table(GWAS_snps).key_by('locus','alleles').repartition(500, shuffle=False)
#    ukb_snps = ukb_snps.annotate(sort_al=hl.sorted(ukb_snps.alleles))
    if verbose:
        print("\nCount 1: " + str(ukb_snps.count()) + '\n')
    
    ukb_snps = ukb_snps.filter(
                    hl.is_snp(ukb_snps.alleles[0], ukb_snps.alleles[1]) &
                    (~(ukb_snps.locus.contig == 'X')) &
                    (~( (ukb_snps.locus.contig == '6') & (ukb_snps.locus.position > 25000000) & (ukb_snps.locus.position < 34000000)))
                )
    if verbose:
        print("\nCount 2: " + str(ukb_snps.count()) + '\n')
    
    # merge in, filter on MAF from the UKBB GWAS sample
    ukb_qc = hl.import_table(GWAS_qc)
    ukb_qc = ukb_qc.annotate(vstruct = hl.parse_variant(ukb_qc.variant))
    ukb_qc = ukb_qc.annotate(locus = ukb_qc.vstruct.locus, alleles = ukb_qc.vstruct.alleles).key_by('locus','alleles')
    ukb_qc2 = ukb_snps.join(ukb_qc.select(ukb_qc.minor_AF))
    if verbose:
        print("\nCount 3: " + str(ukb_qc2.count()) + '\n')
    
    ukb_qc2 = ukb_qc2.filter(
                    (hl.float(ukb_qc2.minor_AF) > 0.01) & 
                    (hl.float(ukb_qc2.minor_AF) < 0.99)
                )
    if verbose:
        print("\nCount 4: " + str(ukb_qc2.count()) + '\n')
    
    # merge in rsid, info (from full UKB sample)
    # and filter to info > 0.9
    ukb_mfi = hl.read_table(GWAS_mfi).key_by('locus','alleles').repartition(500, shuffle=False)
    ukb_qc3 = ukb_qc2.join(ukb_mfi.select('varid','variant','rsid','info'))
    if verbose:
        print("\nCount 5: " + str(ukb_qc3.count()) + '\n')
    
    ukb_qc3 = ukb_qc3.filter(
                    (ukb_qc3.info > 0.9)
                )
    if verbose:
        print("\nCount 6: " + str(ukb_qc3.count()) + '\n')
    
    # drop multi-allelic sites
    loc_count = (ukb_qc3.group_by(ukb_qc3.locus)
                        .aggregate(nloc=hl.agg.count())
                )
    loc_count = loc_count.filter(loc_count.nloc==1)
    if verbose:
        print("\nCount 7: " + str(loc_count.count()) + '\n')
    
    ukb_qc3 = ukb_qc3.key_by('locus').join(loc_count).drop('nloc')
    if verbose:
        print("\nCount 8: " + str(ukb_qc3.count()) + '\n')
    
    # get full hm3 list (autosomal, biallelic snps)
    # and drop extra info
    hm3_snps = hl.read_table(BUCKET + HM3_bi_af).key_by('locus','alleles','rsid').repartition(500, shuffle=False)
    hm3_snps = hm3_snps.select()
    if verbose:
        print("\nCount 9: " + str(hm3_snps.count()) + '\n')

    # merge UKB SNPs to HM3 sites
    # - matching on locus + alleles + rsid
    #    -- are 1-2k SNPs with matching locus and different rsid
    #.   -- trusting ukb and hm3 ref/alt order to match
    # - filter strand ambiguous 
    ukb_hm3 = hm3_snps.join(ukb_qc3.key_by('locus','alleles','rsid'))
    if verbose:
        print("\nCount 10: " + str(ukb_hm3.count()) + '\n')
    
    ukb_hm3 = ukb_hm3.filter(
                        (~hl.is_strand_ambiguous(ukb_hm3.alleles[0], ukb_hm3.alleles[1]))
                    )
    if verbose:
        print("\nCount 11: " + str(ukb_hm3.count()) + '\n')
    
#    hm3_snps = hm3_snps.annotate(sort_al=hl.sorted(hm3_snps.alleles)).key_by('locus','sort_al','rsid')
#    ukb_snps = ukb_snps.annotate(sort_al=hl.sorted(ukb_snps.alleles)).key_by('locus','sort_al','rsid')
#    ukb_hm3 = ukb_snps.join(hm3_snps).key_by('locus','alleles','rsid')
#    ukb_hm3 = ukb_hm3.drop(ukb_hm3.alleles_1)
#    print("\nCount 10: " + str(ukb_hm3.count()) + '\n')
    
    # save passing
#    ukb_qc2 = ukb_qc2.drop(ukb_qc2.sort_al)
    ukb_hm3.write(BUCKET + HM3_qcpos_stem + '.ht', overwrite=True)



if gs_missing(BUCKET + HM3_qcpos_stem + '.tsv.bgz') or force_all:

    # format tsv version
    ukb_hm3 = hl.read_table(BUCKET + HM3_qcpos_stem + '.ht')
    ukb_hm3 = ukb_hm3.annotate(A1=ukb_hm3.alleles[0],A2=ukb_hm3.alleles[1])
    ukb_hm3 = ukb_hm3.key_by('variant')
    ukb_hm3 = ukb_hm3.select(ukb_hm3.locus,ukb_hm3.A1,ukb_hm3.A2,ukb_hm3.rsid,ukb_hm3.varid,UKB_GWAS_MAF=ukb_hm3.minor_AF,UKB_all_INFO=ukb_hm3.info)
    
    # save
    ukb_hm3.export(BUCKET + HM3_qcpos_stem + '.tsv.bgz')



####
#
# Check success and exit
#
####
    
if not gs_missing(BUCKET + HM3_qcpos_stem + '.ht'):

    end_message='''
    ##############
    Finished!!
    
    HapMap3 source:
    {}
    
    HapMap3 raw vcf:
    {}
    
    HapMap3 as Hail MatrixTable: 
    {}
    
    HapMap3 autosomal, biallelic SNPs with freqs table:
    {}
    
    HapMap3 SNPs passing UKBB GWAS QC:
    {}
    {}
    ##############
    '''.format( HM3_FTP + HM3_NAME,
                BUCKET + HM3_NAME,
                BUCKET + HM3_vds,
                BUCKET + HM3_bi_af,
                BUCKET + HM3_qcpos_stem + '.ht',
                BUCKET + HM3_qcpos_stem + '.tsv.bgz'
            )
    
    print(end_message)

else:
    
    print("!!! Error: Failed to write? !!!")

# eof
