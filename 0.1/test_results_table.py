
import hail as hl

ht = hl.read_table('gs://ukb31063-mega-gwas/results-0.2-tables/gwas.imputed_v3.female.ht')
ht = ht.filter(~hl.is_nan(ht.beta))
ht.describe()
ht.show()

