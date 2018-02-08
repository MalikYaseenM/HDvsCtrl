#import pandas as pd


# phenotype_data = 'phs000424.v7.pht002742.v7.p2.c1.GTEx_Subject_Phenotypes.GRU.txt'
# meta_file = pd.read_csv('sample_info_design.csv', sep=',')



# #GTEX-1117F
# #[0:9]
# #GTEX-11GSP_BA9_mRNASeq
# df =  pd.read_csv(phenotype_data, sep='\t', comment='#')
# df = df.head()
# meta_file = meta_file.head()
# meta_file.to_csv('shortsample.csv', sep=',')
# df.to_csv('shortpheno.tsv', sep='\t')

# df = df[['SEX','AGE','SUBJID']]
# SUBJID = df['SUBJID'].values.tolist()
# df = df.set_index([SUBJID])
# df = df.drop(['SUBJID'], axis=1)
# #print(df.head())
# print(df.loc(['GTEX-1117F'],).tolist())




#print(df.head())
#print(meta_file.head())
#print(df.head())
# print(meta_file['identifier'])
#df = pd.read_csv('../../samples/GTEx/all_salmon_quant.tsv')




#print(df['SUBJID'].head())
