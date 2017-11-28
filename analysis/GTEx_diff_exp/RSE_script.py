import os
import pandas as pd

loc = "../../samples/GTEx/results/RSE_firth.tsv"
df = pd.read_csv(loc, sep='\t')

counts_p = df[df['counts__padj'] < .001]
print(counts_p)
