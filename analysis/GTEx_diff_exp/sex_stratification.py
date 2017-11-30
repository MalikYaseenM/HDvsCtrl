import pandas as pd, os
#adds gene id to the usecols list
foo = ['gene_id']

# list with gene_id and female or male subjects
male = foo + pd.read_csv(os.path.abspath("../../samples/GTEx/results/male_sample_design_info.csv"),usecols=['identifier'])['identifier'].tolist()
female = foo + pd.read_csv(os.path.abspath("../../samples/GTEx/results/female_sample_design_info.csv"),usecols=['identifier'])['identifier'].tolist()

# read in the normalized counts use only the males
male_df = pd.read_csv(os.path.abspath("../../samples/GTEx/results/GTEX_salmon_norm.csv"), sep=',', usecols=male)

# read in the normalized counts use only the females
female_df = pd.read_csv(os.path.abspath("../../samples/GTEx/results/GTEX_salmon_norm.csv"), sep=',', usecols=female)

# write out seperate files for males and females
male_df.to_csv(os.path.abspath("../../samples/GTEx/results/GTEX_male_salmon_norm.csv"), sep=',', index=False)
female_df.to_csv(os.path.abspath("../../samples/GTEx/results/GTEX_female_salmon_norm.csv"), sep=',', index=False)
