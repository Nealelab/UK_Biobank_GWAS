###########
# Creates hail table with UKB-provided mfi data (maf, info) for GWAS variants
###########

import hail as hl
hl.init()


######
# Load MFI files
######

# init chr 1
# note: can't wildcard load because of needing to annotate in contig names
ht_mfi = hl.import_table('gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_mfi_chr1_v3.txt', no_header=True)
ht_mfi = ht_mfi.annotate(chrom = '1')

# loop through autosomes
for i in range(2,23):
    ht_mfi_tmp = hl.import_table('gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_mfi_chr'+str(i)+'_v3.txt', no_header=True)
    ht_mfi_tmp = ht_mfi_tmp.annotate(chrom = str(i))
    ht_mfi = ht_mfi.union(ht_mfi_tmp)

# get chrX and PAR
ht_mfi_tmp = hl.import_table('gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_mfi_chrX_v3.txt', no_header=True)
ht_mfi_tmp = ht_mfi_tmp.annotate(chrom = 'X')
ht_mfi = ht_mfi.union(ht_mfi_tmp)

ht_mfi_tmp = hl.import_table('gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_mfi_chrXY_v3.txt', no_header=True)
ht_mfi_tmp = ht_mfi_tmp.annotate(chrom = 'X')
ht_mfi = ht_mfi.union(ht_mfi_tmp)

# add column names
ht_mfi = ht_mfi.rename({'f0': 'varid',
                        'f1': 'rsid',
                        'f2': 'position',
                        'f3': 'allele1_ref',
                        'f4': 'allele2_alt',
                        'f5': 'maf',
                        'f6': 'minor_allele',
                        'f7': 'info'})


######
# set up variant ID for merge with GWAS variant list
######

# construct properly formatted variant name
ht_mfi = ht_mfi.annotate(variant = hl.delimit(hl.array([
                                                        ht_mfi['chrom'],
                                                        hl.str(ht_mfi['position']),
                                                        ht_mfi['allele1_ref'],
                                                        ht_mfi['allele2_alt']]),
                                              delimiter=':'))

# prep to merge with GWAS variant list
ht_mfi = ht_mfi.key_by('variant')
ht_mfi = ht_mfi.annotate(maf = hl.float(ht_mfi.maf), info = hl.float(ht_mfi.info))
ht_mfi = ht_mfi.select('varid', 'rsid', 'maf', 'info')


#######
# load GWAS variant list
#######

# get GWAS variant list
ht_sites = hl.read_table('gs://ukb31063/ukb31063.neale_gwas_variants.ht')
ht_sites = ht_sites.annotate(variant = hl.variant_str(ht_sites.locus, ht_sites.alleles))
ht_sites = ht_sites.key_by('variant')


########
# merge and save
########

# get final merged file with maf/info of the gwas variants
ht = ht_mfi.join(ht_sites, how='inner')
ht = ht.select('locus', 'alleles', 'varid', 'rsid', 'maf', 'info')
print(ht.count())

# save both ht and tsv
ht.write('gs://ukb31063/ukb31063.neale_gwas_variants.imputed_v3.mfi.ht', overwrite=True)
ht.export('gs://ukb31063/ukb31063.neale_gwas_variants.imputed_v3.mfi.tsv.bgz')

###
# for debug: find gwas sites failing to merge
# ht = ht_sites.anti_join(ht_mfi)
# ht.count()
# ht.export('gs://ukb31063/ukb31063.neale_gwas_variants.imputed_v3.missing_mfi.tsv.bgz')
###
