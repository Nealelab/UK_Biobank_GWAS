
from hail import *
hc = HailContext()

kt_mfi = hc.import_table('gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_mfi_chr*_v3.txt', no_header=True)
kt_mfi = kt_mfi.rename({'f0': 'varid',
                        'f1': 'rsid',
                        'f2': 'position',
                        'f3': 'allele1_ref',
                        'f4': 'allele2_alt',
                        'f5': 'maf',
                        'f6': 'minor_allele',
                        'f7': 'info'})
kt_mfi = kt_mfi.key_by('varid')
kt_mfi = kt_mfi.annotate('mfi = {maf: maf.toFloat(), info: info.toFloat()}')
kt_mfi = kt_mfi.select(['varid', 'mfi'])

kt_sites = hc.read('gs://ukb31063-mega-gwas/qc/ukb31063.imputed_v3.sites.vds').variants_table().flatten()
kt_sites = kt_sites.rename({'va.varid': 'varid', 'va.rsid': 'rsid'})
kt_sites = kt_sites.select(['varid', 'rsid', 'v'])
kt_sites = kt_sites.key_by('varid')

kt = kt_mfi.join(kt_sites, how='inner')
kt = kt.key_by('v')
vds = VariantDataset.from_table(kt)
vds.write('gs://ukb31063-mega-gwas/qc/ukb31063.imputed_v3.mfi.vds', overwrite=True)
