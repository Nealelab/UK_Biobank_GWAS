
import hail as hl

ht = hl.import_table('gs://phenotype_31063/ukb31063.sample_qc.tsv.bgz',
                     key='iid')