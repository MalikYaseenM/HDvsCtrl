import pandas as pd
import os

fn = os.path.abspath('../HD_mRNASeq_sample_info.csv')
files = pd.read_csv(fn, sep=',')
print(files)
