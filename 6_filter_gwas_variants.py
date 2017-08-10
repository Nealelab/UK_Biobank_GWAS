from __future__ import print_function
from pprint import pprint
from hail import *

hc = HailContext()

BUCKET = 'gs://ukbb_association/'
APPLICATION = 'test'

VARIANT_VDS = BUCKET + 'all_variants.vds'

GWAS_VARIANTS_VDS = BUCKET + 'gwas_variants.vds'
GWAS_VARIANTS_TSV = BUCKET + 'gwas_variants.tsv'

vds = (
    hc
    .read(VARIANT_VDS)
    .filter_variants_expr("""
        va.isHRC &&
        va.info > 0.8 &&
        va.qc.AF > 0.001 &&
        va.qc.AF < 0.999 &&
        va.qc.pHWE > 1e-10 && 
        va.qc.callRate > 0.95
    """)
)

vds.write(GWAS_VARIANTS_VDS, overwrite=True)

vds = hc.read(GWAS_VARIANTS_VDS)
n_variants = vds.count_variants()

print('')
print('nVariants: ', '{:,}'.format(n_variants))
pprint(vds.variant_schema)

(
    vds
    .export_variants(
        output=GWAS_VARIANTS_TSV,
        expr="""
            variant = Variant(v.contig.replace("^0",""), v.start, v.ref, v.alt()),
            chr = v.contig.replace("^0",""),
            pos = v.start,
            ref = v.ref,
            alt = v.alt(),
            rsid = va.rsid,
            info = va.info,
            va.qc.*
        """
    )
)
