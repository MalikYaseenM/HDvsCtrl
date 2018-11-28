import random, pandas as pd
from os import path
sample_info = path.join(path.abspath('../../samples/GTEx/results/'),
                        'new_sample_info.csv')
path_name = path.abspath('../../samples/GTEx/results/sample_info')

df = pd.read_csv(sample_info, sep=',')
for idx in range(1000):
    new_df = df
    random.shuffle(new_df['brain_region'])    
    new_df_name=path_name + '_V' + str(idx) +'.csv'
    new_df.to_csv(new_df_name, sep=',', index=False)
    
    
# new_df = df
# random.shuffle(new_df['brain_region'])
# print(new_df)



