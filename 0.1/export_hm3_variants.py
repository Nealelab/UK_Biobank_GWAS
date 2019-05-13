
from hail import *
hc = HailContext()

kt = hc.read_table('gs://ukb31063-mega-gwas/ldsc/ld_ref_panel/hm3.r3.b37.auto_bi_af.ukbb_gwas_qcpos.kt')
kt.export('gs://ukb31063-mega-gwas/hm3_variants.tsv.bgz')
