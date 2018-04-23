
from hail import *
hc = HailContext()

vds = hc.read('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.imputed_v3.sites.vds').repartition(2000)
vds = vds.vep(config='/vep/vep-gcloud.properties')
vds.write('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.imputed_v3.vep.vds', overwrite=True)
