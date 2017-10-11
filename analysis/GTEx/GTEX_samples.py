import csv
import pandas
#usecols=col_names)
col_names = ['identifier', 'sample_name', 'brain_region', 'abspath']
df = pandas.read_csv('GTEX_info.tsv', sep='\t', usecols=col_names)

samples = df['identifier'].tolist()
