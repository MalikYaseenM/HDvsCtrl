import pandas as pd

## List of bad genes so far:
# Gene name, Filtered file, Line from all_salmon_norm_filter
# ENSG00000280594.1, Filtered, line 237

df = pd.read_csv("all_salmon_norm_filter.csv")

# drops row 0 to 236 returns rows starting from 237
# drop bad gene #1
df = df.drop(df.index[0:236])




# write to csv file
df.to_csv("clean_filtered_genes.csv", index=False)

