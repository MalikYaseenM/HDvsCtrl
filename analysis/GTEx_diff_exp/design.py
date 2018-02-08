import os, pandas as pd
# CODE 1 is male CODE 2 is female
fn = os.path.abspath('../GTEx/GTEX_info.tsv')

info_file = pd.DataFrame(pd.read_csv(fn, sep='\t',usecols=['identifier','brain_region']))

info_file.to_csv(os.path.abspath('../../samples/GTEx/results/sample_info_design.csv'), sep=',',index=False)
