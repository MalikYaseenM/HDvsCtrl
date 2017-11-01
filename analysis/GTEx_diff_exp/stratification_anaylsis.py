import pandas as pd

fem_df = pd.read_csv('GTEX_female_salmon_firth.tsv', sep='\t')
m_df = pd.read_csv('GTEX_male_salmon_firth.tsv', sep='\t')


# print(fem_df[fem_df['counts__padj'] < .05].shape)
print(m_df[m_df['counts__padj'] < .05].shape)
