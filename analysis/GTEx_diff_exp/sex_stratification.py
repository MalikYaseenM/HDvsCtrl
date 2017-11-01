import pandas as pd
foo = ['gene_id']
male = foo + pd.read_csv('male_sample_design_info.csv',
                         usecols=['identifier'])['identifier'].tolist()
female = foo + pd.read_csv('female_sample_design_info.csv',
                           usecols=['identifier'])['identifier'].tolist()

male_df = pd.read_csv('GTEX_salmon_norm.csv', sep=',', usecols=male)
female_df = pd.read_csv('GTEX_salmon_norm.csv', sep=',', usecols=female)

male_df.to_csv("GTEX_male_salmon_norm.csv", sep=',', index=False)
female_df.to_csv("GTEX_female_salmon_norm.csv", sep=',', index=False)
