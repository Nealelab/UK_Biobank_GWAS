
from hail import *
hc = HailContext()

(hc.read('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.imputed_v3.mfi.vds')
   .variants_table()
   .filter('va.mfi.info >= 0.8 && v.contig == "X"', keep=True)
   .select('v')
   .repartition(50)
   .write('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.mfi.chrX.filterRepartition.kt', overwrite=True))

(hc.read('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.imputed_v3.variant_qc.chrX.vds')
   .variants_table()
   .filter("""va.qc.AF >= 0.001 &&
              va.qc.AF <= 0.999 &&
              va.qc.pHWE > 1e-10 &&
              va.qc.callRate >= 0.95""", keep=True)
   .select('v')
   .repartition(50)
   .write('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.qc.chrX.filterRepartition.kt', overwrite=True))

(hc.read('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.imputed_v3.consequences_categorized.vep.vds')
   .variants_table()
   .filter('["ptv", "missense", "synonymous"].toSet().contains(va.consequence_category) && v.contig == "X"', keep=True)
   .select('v')
   .repartition(50)
   .write('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.vep.chrX.filterRepartition.kt', overwrite=True))

kt_mfi = hc.read_table('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.mfi.chrX.filterRepartition.kt')
kt_qc = hc.read_table('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.qc.chrX.filterRepartition.kt')
kt_vep = hc.read_table('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.vep.chrX.filterRepartition.kt')

kt_qc_info = kt_qc.join(kt_mfi, how='inner')
kt_vep_info = kt_vep.join(kt_mfi, how='inner')
kt_join = kt_qc_info.join(kt_vep_info, how='outer')

vds = VariantDataset.from_table(kt_join).annotate_variants_expr('va.keep = true').cache()

(hc.read('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.imputed_v3.mfi.vds')
   .annotate_variants_vds(vds, expr='va.keep = vds.keep')
   .filter_variants_expr('isDefined(va.keep)')
   .repartition(50)
   .write('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.mfi.gwas_variants.chrX.vds', overwrite=True))
vds_mfi = hc.read('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.mfi.gwas_variants.chrX.vds')

(hc.read('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.imputed_v3.variant_qc.chrX.vds')
   .annotate_variants_vds(vds, expr='va.keep = vds.keep')
   .filter_variants_expr('isDefined(va.keep)')
   .repartition(50)
   .write('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.qc.gwas_variants.chrX.vds', overwrite=True))
vds_qc = hc.read('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.qc.gwas_variants.chrX.vds')

(hc.read('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.imputed_v3.consequences_categorized.vep.vds')
   .annotate_variants_vds(vds, expr='va.keep = vds.keep')
   .filter_variants_expr('isDefined(va.keep)')
   .repartition(50)
   .write('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.vep.gwas_variants.chrX.vds', overwrite=True))
vds_vep = hc.read('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.vep.gwas_variants.chrX.vds')

vds = (vds.annotate_variants_vds(vds_mfi, expr='va.rsid = vds.rsid, va.varid = vds.varid, va.info = vds.mfi.info')
          .annotate_variants_vds(vds_qc, expr='va.qc = vds.qc')
          .annotate_variants_vds(vds_vep, root='va.vep')
          .filter_variants_expr('isDefined(va.qc.AF)'))

(vds.variants_table()
    .select('v')
    .export('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.gwas_variants.chrX.tsv'))

vds.write('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.gwas_variants.chrX.vds', overwrite=True)
print 'nVariants: {:,}'.format(hc.read('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.gwas_variants.chrX.vds').count_variants())
