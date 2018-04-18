
from hail import *
hc = HailContext()

(hc.read('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.imputed_v3.mfi.vds')
   .variants_table()
   .filter('va.mfi.info >= 0.8 && v.contig != "X"', keep=True)
   .annotate('v = let c = if (v.contig.length() == 1) "0" + v.contig else v.contig in Variant(c, v.start, v.ref, v.alt())')
   .select('v')
   .repartition(1000)
   .write('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.mfi.autosomes.filterRepartition.kt', overwrite=True))

(hc.read('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.imputed_v3.variant_qc.autosomes.vds')
   .variants_table()
   .filter("""va.qc.AF >= 0.001 &&
              va.qc.AF <= 0.999 &&
              va.qc.pHWE > 1e-10 &&
              va.qc.callRate >= 0.95""", keep=True)
   .select('v')
   .repartition(1000)
   .write('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.qc.autosomes.filterRepartition.kt', overwrite=True))

(hc.read('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.imputed_v3.consequences_categorized.vep.vds')
   .variants_table()
   .filter('["ptv", "missense", "synonymous"].toSet().contains(va.consequence_category) && v.contig != "X"', keep=True)
   .annotate('v = let c = if (v.contig.length() == 1) "0" + v.contig else v.contig in Variant(c, v.start, v.ref, v.alt())')
   .select('v')
   .repartition(1000)
   .write('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.vep.autosomes.filterRepartition.kt', overwrite=True))

kt_mfi = hc.read_table('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.mfi.autosomes.filterRepartition.kt')
kt_qc = hc.read_table('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.qc.autosomes.filterRepartition.kt')
kt_vep = hc.read_table('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.vep.autosomes.filterRepartition.kt')

kt_qc_info = kt_qc.join(kt_mfi, how='inner')
kt_vep_info = kt_vep.join(kt_mfi, how='inner')
kt_join = kt_qc_info.join(kt_vep_info, how='outer')

VariantDataset.from_table(kt_join).write('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.gwas_variants.autosomes.vds', overwrite=True)

vds = hc.read('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.gwas_variants.autosomes.vds').annotate_variants_expr('va.keep = true')
print 'nVariants: {:,}'.format(vds.count_variants())

(vds.variants_table()
    .select('v')
    .annotate('v = Variant(v.contig.replace("^0", ""), v.start, v.ref, v.alt())')
    .export('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.gwas_variants.autosomes.tsv'))

(VariantDataset.from_table(hc.read('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.imputed_v3.mfi.vds')
                             .variants_table()
                             .annotate('v = let c = if (v.contig.length() == 1 && v.contig != "X") "0" + v.contig else v.contig in Variant(c, v.start, v.ref, v.alt())'))
               .annotate_variants_expr('va = va.va')
               .annotate_variants_vds(vds, expr='va.keep = vds.keep')
               .filter_variants_expr('isDefined(va.keep)')
               .repartition(1000)
               .write('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.mfi.gwas_variants.autosomes.vds', overwrite=True))
vds_mfi = hc.read('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.mfi.gwas_variants.autosomes.vds')

(hc.read('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.imputed_v3.variant_qc.autosomes.vds')
   .annotate_variants_vds(vds, expr='va.keep = vds.keep')
   .filter_variants_expr('isDefined(va.keep)')
   .repartition(1000)
   .write('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.qc.gwas_variants.autosomes.vds', overwrite=True))
vds_qc = hc.read('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.qc.gwas_variants.autosomes.vds')

(VariantDataset.from_table(hc.read('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.imputed_v3.consequences_categorized.vep.vds')
                             .variants_table()
                             .annotate('v = let c = if (v.contig.length() == 1 && v.contig != "X") "0" + v.contig else v.contig in Variant(c, v.start, v.ref, v.alt())'))
               .annotate_variants_expr('va = va.va')
               .annotate_variants_vds(vds, expr='va.keep = vds.keep')
               .filter_variants_expr('isDefined(va.keep)')
               .repartition(100)
               .write('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.vep.gwas_variants.autosomes.vds', overwrite=True))
vds_vep = hc.read('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.vep.gwas_variants.autosomes.vds')

vds = (vds.annotate_variants_vds(vds_mfi, expr='va.rsid = vds.rsid, va.varid = vds.varid, va.info = vds.mfi.info')
          .annotate_variants_vds(vds_qc, expr='va.qc = vds.qc')
          .annotate_variants_vds(vds_vep, root='va.vep'))

vds.write('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.gwas_variants.autosomes.with_qc_annotations.vds', overwrite=True)
