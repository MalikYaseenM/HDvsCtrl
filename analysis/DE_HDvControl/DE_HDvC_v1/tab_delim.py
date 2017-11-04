import pandas as pd

df = pd.read_csv("clean_salmon_norm_firth.csv", sep="\t")
df.to_csv("clean_salmon_norm_firth_tab.csv", index=False)
