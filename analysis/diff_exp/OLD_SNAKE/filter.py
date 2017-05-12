import pandas as pd

#df = pd.read_csv("head_all_salmon_quant_norm.csv")
df = pd.read_csv("../../samples/all_salmon_quant_norm.csv")

df['means'] = df.mean(axis=1)
df = df[df.means > 10]
df = df.drop('means', axis=1)

#df.to_csv("head_all_salmon_norm_filter.csv", index=False)
df.to_csv("all_salmon_norm_filter.csv", index=False)
