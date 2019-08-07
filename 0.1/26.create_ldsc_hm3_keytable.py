#! /usr/bin/env python

####
# setup
####

# re-run everything?
force_all = False
force_final = True

### input

###
### NOTE: these paths are now deprecated
###

# where to save files
BUCKET = 'gs://ukb31063-mega-gwas/ldsc/ld_ref_panel/'

# source of HapMap3 sites file (path only)
HM3_FTP = 'ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/'

# name of original HapMap3 sites file
HM3_NAME = 'hapmap_3.3_b37_pop_stratified_af.vcf.gz'

# UKBB variant QC in the GWAS samples (for MAF, info, passing list)
GWAS_qc = 'gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.gwas_variants.autosomes.vds'

### output

# name to save HapMap3 sites as matrixtable
HM3_vds = 'hm3.r3.b37.vds'

# name to save autosomal, biallelic SNPs (also has parsed population allele freqs)
HM3_bi_af = 'hm3.r3.b37.auto_bi_af.kt'

# name to save HapMap3 SNPs passing ldsc QC in UKBB GWAS sample
# - will save both Hail table (.kt) and flat file (.tsv)
HM3_qcpos_stem = 'hm3.r3.b37.auto_bi_af.ukbb_gwas_qcpos.no_mhc'



# load packages
print "Preparing packages, etc..."
from hail import *
import subprocess
hc = HailContext()


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

    print "Downloading HapMap3 vcf..."
    
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
    
    print "Creating Hail table of HapMap3 sites..."
        
    (hc
        .import_vcf(BUCKET + HM3_NAME, force_bgz=True, min_partitions=500)
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
    print "Creating keytable of autosomal, biallelic HM3 SNPs with population allele freqs..."
    ###
    
    # read full HM3 sites
    hm3_vds = hc.read(BUCKET + HM3_vds)
    
    # filter to autosomal, biallelic snps and extract allelic freqs
    hm3_keep = (hm3_vds.filter_variants_expr('''
                                                v.isBiallelic() && 
                                                v.isAutosomal() &&
                                                v.altAllele().isSNP()
                                            ''')    
                        .annotate_variants_expr('''va.ASW_AF = 1.0 - va.info.ASW[0].split("=")[1].toFloat(), 
                                                   va.CEU_AF = 1.0 - va.info.CEU[0].split("=")[1].toFloat(), 
                                                   va.CHB_AF = 1.0 - va.info.CHB[0].split("=")[1].toFloat(),
                                                   va.CHD_AF = 1.0 - va.info.CHD[0].split("=")[1].toFloat(),
                                                   va.CHS_AF = 1.0 - va.info.CHS[0].split("=")[1].toFloat(),
                                                   va.CLM_AF = 1.0 - va.info.CLM[0].split("=")[1].toFloat(),
                                                   va.FIN_AF = 1.0 - va.info.FIN[0].split("=")[1].toFloat(),
                                                   va.GBR_AF = 1.0 - va.info.GBR[0].split("=")[1].toFloat(),
                                                   va.GIH_AF = 1.0 - va.info.GIH[0].split("=")[1].toFloat(),
                                                   va.IBS_AF = 1.0 - va.info.IBS[0].split("=")[1].toFloat(),
                                                   va.JPT_AF = 1.0 - va.info.JPT[0].split("=")[1].toFloat(),
                                                   va.LWK_AF = 1.0 - va.info.LWK[0].split("=")[1].toFloat(),
                                                   va.MKK_AF = 1.0 - va.info.MKK[0].split("=")[1].toFloat(),
                                                   va.MXL_AF = 1.0 - va.info.MXL[0].split("=")[1].toFloat(),
                                                   va.PUR_AF = 1.0 - va.info.PUR[0].split("=")[1].toFloat(),
                                                   va.TSI_AF = 1.0 - va.info.TSI[0].split("=")[1].toFloat(),
                                                   va.YRI_AF = 1.0 - va.info.YRI[0].split("=")[1].toFloat()
                                                '''))


    # save keytable
    hm3_keep.variants_table().write(BUCKET + HM3_bi_af, overwrite=True)


###
# Create filtered keytable with autosomal, biallelic HM3 snps
# with MAF > .01, INFO > 0.9 and passing QC in UKBB GWAS analysis set
#
# Also provides both UKBB and HM3 formatting of variant
###

if gs_missing(BUCKET + HM3_qcpos_stem + '.kt') or force_all or force_final:
    
    print "Creating Hail table of HM3 SNPs passing UKBB GWAS QC..."


    # get list of SNPs to be used in GWAS with QC info
    # apply MAF, info, SNP filters here pre-merge for efficiency
    ukb_snps_all = hc.read(GWAS_qc)
    print ukb_snps_all.count()
    ukb_snps = (ukb_snps_all.variants_table()
                            .flatten()
                            .filter('''
                                        `va.keep` &&
                                        `va.info` > 0.9 &&
                                        `va.qc.AF` > .01 &&
                                        `va.qc.AF` < .99 &&
                                        v.isBiallelic() &&
                                        v.altAllele().isSNP() &&
                                        !(v.contig == "06" && v.start > 25000000 && v.start < 34000000)
                                    ''')
                            .annotate('''v_hm3 = Variant(v.contig.replace("^0",""), v.start, v.ref, v.alt), 
                                         v_ukb = v,
                                         locus = Locus(v.contig.replace("^0",""), v.start)''')
                            .select(['v_ukb','v_hm3','locus','va.keep','va.rsid','va.varid','va.info','va.qc.AF'])
                            
                )
    print ukb_snps.count()
    
    # drop multi-allelic loci
    # - required to avoid issues with allele discordance 
    #   between HM3 and the w_hm3.snplist
    loc_count = (ukb_snps.key_by('locus')
                         .aggregate_by_key('locus = locus', 'num_var = v_hm3.count().toInt()')
                         .filter('num_var == 1')
                )
    ukb_snps_bi = ukb_snps.key_by('locus').join(loc_count).drop(['num_var'])
    print ukb_snps_bi.count()
    
    # get full hm3 list (autosomal, biallelic snps)
    # and drop extra info
    hm3_snps_all = hc.read_table(BUCKET + HM3_bi_af)
    hm3_snps = hm3_snps_all.flatten().select(['v','va.rsid']).rename({'va.rsid': 'hm3_rsid'})
    
    # merge UKB SNPs to HM3 sites on variant (chr:pos:ref:alt)
    merge_snps = hm3_snps.key_by('v').join(ukb_snps_bi.key_by('v_hm3'))
    print merge_snps.count()
    
    # filter to rsid match, not strand ambi
    keep_snps = (merge_snps.filter('''
                                     `va.rsid` == hm3_rsid &&
                                     !(v.ref=='A' && v.alt=='T') &&
                                     !(v.ref=='T' && v.alt=='A') &&
                                     !(v.ref=='C' && v.alt=='G') &&
                                     !(v.ref=='G' && v.alt=='C')                                    
                                ''')
                            .drop(['va.keep','hm3_rsid','va.varid'])
                            .rename({
                              'va.rsid':'rsid',
                              'va.info':'UKBB_INFO',
                              'va.qc.AF':'UKBB_GWAS_AF'
                            })
                )
    print keep_snps.count()

    # save passing
    keep_snps.write(BUCKET + HM3_qcpos_stem + '.kt', overwrite=True)



if gs_missing(BUCKET + HM3_qcpos_stem + '.tsv.bgz') or force_all or force_final:

    # format tsv version
    keep_snps = hc.read_table(BUCKET + HM3_qcpos_stem + '.kt')
    
    # save
    keep_snps.export(BUCKET + HM3_qcpos_stem + '.tsv.bgz')



####
#
# Check success and exit
#
####
    
if not gs_missing(BUCKET + HM3_qcpos_stem + '.kt'):

    end_message='''
    ##############
    Finished!!
    
    HapMap3 source:
    {}
    
    HapMap3 raw vcf:
    {}
    
    HapMap3 as Hail MatrixTable: 
    {}
    
    HapMap3 autosomal, biallelic SNPs with freqs keytable:
    {}
    
    HapMap3 SNPs passing UKBB GWAS QC:
    {}
    {}
    ##############
    '''.format( HM3_FTP + HM3_NAME,
                BUCKET + HM3_NAME,
                BUCKET + HM3_vds,
                BUCKET + HM3_bi_af,
                BUCKET + HM3_qcpos_stem + '.kt',
                BUCKET + HM3_qcpos_stem + '.tsv.bgz'
            )
    
    print end_message

else:
    
    print "!!! Error: Failed to write? !!!" 

# eof
