import pandas as pd

# Creates a DataFrame which contains GTEX id and SRR id
# Brain - Frontal Cortex (BA9)
# Brain - Caudate (basal ganglia)
# cols = ['Run_s','dataset_id', 'SMTSD']
# df = pd.read_csv('../reference/all_info.csv', sep=',', usecols=cols)
# df = df.loc[df['SMTSD'].isin(['Brain - Frontal Cortex (BA9)', 'Brain - Caudate (basal ganglia)'])]
# sub_df = df[['dataset_id', 'Run_s']]
# sub_df.to_csv('srr&id.tsv', sep='\t', index=False))

# Creates new GTEX_info file which contains columns ['identifier','Run_s', 'sample_name','brain_region','abspath']
# key = pd.read_csv('srr&id.tsv', sep='\t')
# df = pd.read_csv('../GTEx/GTEX_info.tsv', sep='\t')
# ids = df.identifier.tolist()
# sub_df = key[key['dataset_id'].isin(ids)]
# merge_df = df.merge(sub_df, left_on='identifier', right_on='dataset_id', how='outer')
# ultimate_df = merge_df[['identifier','Run_s', 'sample_name','brain_region','abspath']]
# ultimate_df.to_csv('new_GTEX_info.tsv', sep='\t', index=False)

#Converts SRR IDs to GTEX ids 
sample_info = pd.read_csv('new_GTEX_info.tsv', sep='\t')
keys = sample_info['Run_s'].values.tolist()
values = sample_info['identifier'].values.tolist()
key_value = list(zip(keys, values))
sample_dict = dict((x, y) for x, y in key_value)
cols = ['gene_id'] + keys
RSE = pd.read_csv('RSE_Counts.csv', sep=',', usecols=cols)
new_df = RSE.rename(columns=sample_dict)
new_df.to_csv('RSE_ID_Counts.tsv', sep='\t', index=False)
