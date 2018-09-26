import pandas as pd, os

# Merge deseq2 and raw salmon counts
fm = os.path.abspath('../../samples/Analysis_Results/mart_export.txt')
fq = os.path.abspath('../../samples/all_salmon_quant.tsv')
fr = os.path.abspath('../../samples/Analysis_Results/sympBA9_deseq2_filter_halfzero.csv')
fl = os.path.abspath('../HD_mRNASeq_sample_info.csv')

# Dataframes
df = pd.read_csv(fq, sep='\t', comment='#')
deseq2 = pd.read_csv(fr, sep=',', comment='#')
ensemble = pd.read_csv(fm, sep='\t', comment='#')
ensemble = ensemble.rename(columns = {'Gene stable ID':'gene_id'})
sample_info = pd.read_csv(fl, sep=',', comment='#')

# Get sample IDs 
sample_info = sample_info[(sample_info['Dataset.protocol'] != 'TruSeq Ribo-depleted') & (sample_info['Dataset.dataset_id'] != 'H_0014_BA9_mRNASeq')]
samples = sample_info.iloc[:,[14,1]].copy()
samples.columns = ["Data_id","Subject_type"]

hdpos = samples['Data_id'][samples['Subject_type'].str.contains("HDpos")].tolist()
hdpos_BA9 = [ _ for _ in hdpos if "BA9" in _]
hdpos_CAP = [ _ for _ in hdpos if "CAP" in _]
CAP = [ _ for _ in samples['Data_id'] if "CAP" in _]
C_cap = [ _ for _ in CAP if _.startswith('C')]
HD_cap = [ _ for _ in CAP if _.startswith('H')]
symp_BA9 = samples[samples['Data_id'].str.contains('BA9') & (samples['Subject_type'] != 'HDpos')]
C_BA9 = [ _ for _ in symp_BA9['Data_id'] if _.startswith('C')]
HD_BA9 = [ _ for _ in symp_BA9['Data_id'] if _.startswith('H')]

# Process deseq2, drop irrelevant columns
deseq2 = deseq2[['gene_id', 'padj']]

# Merge
merged = pd.merge(deseq2, df, on='gene_id')
merged['gene_id'] = merged['gene_id'].str.split('.').str[0]
merged = pd.merge(ensemble, merged, on='gene_id')
merged = merged.sort_values(['padj'], ascending=[True])

s = merged.drop(['gene_id', 'HGNC symbol', 'padj'], 1)
s = s[list(s)[:1] + C_BA9 + HD_BA9 + hdpos_BA9 + hdpos_CAP + C_cap]
s.to_csv(os.path.abspath("../../samples/Analysis_Results/sympBA9_counts_deseq2.csv"), index=False)
