
from hail import *
hc = HailContext()

(hc.read('gs://ukb31063-mega-gwas/qc/ukb31063.imputed_v3.mfi.vds')
   .variants_table()
   .filter('va.mfi.info >= 0.8 && v.contig == "X"', keep=True)
   .select('v')
   .repartition(12)
   .write('gs://ukb31063-mega-gwas/qc/ukb31063.mfi.chrXY.filterRepartition.kt', overwrite=True))

(hc.read('gs://ukb31063-mega-gwas/qc/ukb31063.imputed_v3.variant_qc.chrXY.vds')
   .variants_table()
   .filter("""va.qc.AF >= 0.001 &&
              va.qc.AF <= 0.999 &&
              va.qc.pHWE > 1e-10 &&
              va.qc.callRate >= 0.95""", keep=True)
   .select('v')
   .repartition(12)
   .write('gs://ukb31063-mega-gwas/qc/ukb31063.qc.chrXY.filterRepartition.kt', overwrite=True))

(hc.read('gs://ukb31063-mega-gwas/qc/ukb31063.imputed_v3.consequences_categorized.vep.vds')
   .variants_table()
   .filter('["ptv", "missense", "synonymous"].toSet().contains(va.consequence_category) && v.contig == "X"', keep=True)
   .select('v')
   .repartition(12)
   .write('gs://ukb31063-mega-gwas/qc/ukb31063.vep.chrXY.filterRepartition.kt', overwrite=True))

kt_mfi = hc.read_table('gs://ukb31063-mega-gwas/qc/ukb31063.mfi.chrXY.filterRepartition.kt')
kt_qc = hc.read_table('gs://ukb31063-mega-gwas/qc/ukb31063.qc.chrXY.filterRepartition.kt')
kt_vep = hc.read_table('gs://ukb31063-mega-gwas/qc/ukb31063.vep.chrXY.filterRepartition.kt')

print kt_mfi.count()
print kt_qc.count()
print kt_vep.count()

import sys
sys.exit()

kt_qc_info = kt_qc.join(kt_mfi, how='inner')
kt_vep_info = kt_vep.join(kt_mfi, how='inner')
kt_join = kt_qc_info.join(kt_vep_info, how='outer')

kt_qc_info.show()
kt_vep_info.show()
kt_join.show()
print kt_qc_info.count()
print kt_vep_info.count()
print kt_join.count()

import sys
sys.exit()

vds = VariantDataset.from_table(kt_join).annotate_variants_expr('va.keep = true').cache()

(hc.read('gs://ukb31063-mega-gwas/qc/ukb31063.imputed_v3.mfi.vds')
   .annotate_variants_vds(vds, expr='va.keep = vds.keep')
   .filter_variants_expr('isDefined(va.keep)')
   .repartition(50)
   .write('gs://ukb31063-mega-gwas/qc/ukb31063.mfi.gwas_variants.chrXY.vds', overwrite=True))
vds_mfi = hc.read('gs://ukb31063-mega-gwas/qc/ukb31063.mfi.gwas_variants.chrXY.vds')

(hc.read('gs://ukb31063-mega-gwas/qc/ukb31063.imputed_v3.variant_qc.chrXY.vds')
   .annotate_variants_vds(vds, expr='va.keep = vds.keep')
   .filter_variants_expr('isDefined(va.keep)')
   .repartition(50)
   .write('gs://ukb31063-mega-gwas/qc/ukb31063.qc.gwas_variants.chrXY.vds', overwrite=True))
vds_qc = hc.read('gs://ukb31063-mega-gwas/qc/ukb31063.qc.gwas_variants.chrXY.vds')

(hc.read('gs://ukb31063-mega-gwas/qc/ukb31063.imputed_v3.consequences_categorized.vep.vds')
   .annotate_variants_vds(vds, expr='va.keep = vds.keep')
   .filter_variants_expr('isDefined(va.keep)')
   .repartition(50)
   .write('gs://ukb31063-mega-gwas/qc/ukb31063.vep.gwas_variants.chrXY.vds', overwrite=True))
vds_vep = hc.read('gs://ukb31063-mega-gwas/qc/ukb31063.vep.gwas_variants.chrXY.vds')

vds = (vds.annotate_variants_vds(vds_mfi, expr='va.rsid = vds.rsid, va.varid = vds.varid, va.info = vds.mfi.info')
          .annotate_variants_vds(vds_qc, expr='va.qc = vds.qc')
          .annotate_variants_vds(vds_vep, root='va.vep')
          .filter_variants_expr('isDefined(va.qc.AF)'))

(vds.variants_table()
    .select('v')
    .export('gs://ukb31063-mega-gwas/qc/ukb31063.gwas_variants.chrXY.tsv'))

vds.write('gs://ukb31063-mega-gwas/qc/ukb31063.gwas_variants.chrXY.vds', overwrite=True)
print 'nVariants: {:,}'.format(hc.read('gs://ukb31063-mega-gwas/qc/ukb31063.gwas_variants.chrXY.vds').count_variants())
