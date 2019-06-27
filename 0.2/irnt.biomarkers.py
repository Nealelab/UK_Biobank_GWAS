
import sys
import hail as hl
import pandas as pd

sex = sys.argv[1]

if sex not in set(['both_sexes', 'female', 'male']):
    raise ValueError('Invalid argument: "sex" must be one of {"both_sexes", "female", "male"}.')

ht = hl.read_table('gs://ukb31063/ukb31063.biomarkers.ht')
ht_samples = hl.read_table(f'gs://ukb31063/ukb31063.neale_gwas_samples.{sex}.ht')

ht = ht.filter(hl.is_defined(ht_samples[ht.s]))

df = ht.to_pandas()
df.index = df['s']
df = df.drop('s', axis=1)

dfp = df.rank()
dfp = (dfp - 0.5) / (~dfp.isnull()).sum()
dfp.columns = [x + '_prob' for x in dfp.columns]

df.columns = [x + '_raw' for x in df.columns]
df = pd.merge(df, dfp, how='inner', left_index=True, right_index=True)
df.loc[:, 's'] = df.index

ht = hl.Table.from_pandas(df, key='s')
ht = ht.annotate(**{x.replace('_prob', '_irnt'): hl.qnorm(ht[x])
                    for x in ht.row_value if x.endswith('_prob')})
ht = ht.annotate(**{x: hl.or_missing(~hl.is_nan(ht[x]), ht[x]) for x in ht.row_value})

ht = ht.drop(*[x for x in ht.row_value if x.endswith('_prob')])
ht.write(f'gs://ukb31063/ukb31063.biomarkers_gwas.{sex}.ht', overwrite=True)
