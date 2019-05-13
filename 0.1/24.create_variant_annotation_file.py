
from hail import *
hc = HailContext()

vds_autosomes = hc.read('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.gwas_variants.autosomes.with_qc_annotations.vds')
vds_chrX = hc.read('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.gwas_variants.chrX.with_qc_annotations.vds')
vds = VariantDataset.union(vds_autosomes, vds_chrX)

vds = vds.filter_variants_expr('isDefined(va.qc.AF)')
vds = vds.annotate_variants_expr("""va.minor_allele = if (va.qc.AF <= 0.5) v.alt() else v.ref,
                                    va.minor_AF = if (va.qc.AF <= 0.5) va.qc.AF else 1.0 - va.qc.AF""")

vds.export_variants('gs://ukb31063-mega-gwas/hail-0.1/annotations/variants.tsv.gz',
    """
    variant = Variant(v.contig.replace("^0", ""), v.start, v.ref, v.alt()),
    contig = v.contig.replace("^0", ""),
    pos = v.start,
    ref = v.ref,
    alt = v.alt(),
    rsid = va.rsid,
    varid = va.varid,
    consequence = va.vep.consequence,
    consequence_category = va.vep.consequence_category,
    info = va.info,
    call_rate = va.qc.callRate,
    AC = va.qc.AC,
    AF = va.qc.AF,
    minor_allele = va.minor_allele,
    minor_AF = va.minor_AF,
    p_hwe = va.qc.pHWE,
    n_called = va.qc.nCalled,
    n_not_called = va.qc.nNotCalled,
    n_hom_ref = va.qc.nHomRef,
    n_het = va.qc.nHet,
    n_hom_var = va.qc.nHomVar,
    n_non_ref = va.qc.nNonRef,
    r_heterozygosity = va.qc.rHeterozygosity,
    r_het_hom_var = va.qc.rHetHomVar,
    r_expected_het_frequency = va.qc.rExpectedHetFrequency
    """
)
