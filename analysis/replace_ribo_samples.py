# we originally had sequenced these samples as ribo-depleted, but found they
# were difficult to compare to the poly-A samples
# we therefore resequenced the samples with a poly-A library prep and now
# are updating the locations of the files appropriately
import pandas
import os

df = pandas.read_csv('HD_mRNASeq_sample_info.csv',skiprows=1)
df.index = df['Dataset.dataset_id']

# C_0091 wasn't in the database before, so it needs new metadata
fields = ['Subject.subject_id','Subject.subject_type','Subject.external_id',
          'Subject.PMI','Subject.death','Subject.sex']
df.loc['C_0091_BA9_mRNASeq',fields] = ['C_0091','Control','B7804','14.12','89.0','F']

new_samples = [
  ['C_0114_CAP_mRNASeq', '/restricted/projectnb/mlhd/mRNA/2018-1-29_Myers/C_0114_CAP_mRNASeq__R1.fastq.gz', '/restricted/projectnb/mlhd/mRNA/2018-1-29_Myers/C_0114_CAP_mRNASeq__R2.fastq.gz'],
  ['C_0076_BA9_mRNASeq', '/restricted/projectnb/mlhd/mRNA/2018-1-29_Myers/C_0076_BA9_mRNASeq__R1.fastq.gz', '/restricted/projectnb/mlhd/mRNA/2018-1-29_Myers/C_0076_BA9_mRNASeq__R2.fastq.gz'],
  ['H_1106_BA9_mRNASeq', '/restricted/projectnb/mlhd/mRNA/2018-1-29_Myers/H_1106_BA9_mRNASeq__R1.fastq.gz', '/restricted/projectnb/mlhd/mRNA/2018-1-29_Myers/H_1106_BA9_mRNASeq__R2.fastq.gz'],
  ['C_0091_BA9_mRNASeq', '/restricted/projectnb/mlhd/mRNA/2018-1-29_Myers/C_0091_BA9_mRNASeq__R1.fastq.gz', '/restricted/projectnb/mlhd/mRNA/2018-1-29_Myers/C_0091_BA9_mRNASeq__R2.fastq.gz'],
  ['H_1105_CAP_mRNASeq', '/restricted/projectnb/mlhd/mRNA/2018-1-29_Myers/H_1105_CAP_mRNASeq__R1.fastq.gz', '/restricted/projectnb/mlhd/mRNA/2018-1-29_Myers/H_1105_CAP_mRNASeq__R2.fastq.gz'],
  ['C_0074_BA9_mRNASeq', '/restricted/projectnb/mlhd/mRNA/2018-1-29_Myers/C_0074_BA9_mRNASeq__R1.fastq.gz', '/restricted/projectnb/mlhd/mRNA/2018-1-29_Myers/C_0074_BA9_mRNASeq__R2.fastq.gz'],
  ['H_1104_BA9_mRNASeq', '/restricted/projectnb/mlhd/mRNA/2018-1-29_Myers/H_1104_BA9_mRNASeq__R1.fastq.gz', '/restricted/projectnb/mlhd/mRNA/2018-1-29_Myers/H_1104_BA9_mRNASeq__R2.fastq.gz'],
  ['H_1105_BA9_mRNASeq', '/restricted/projectnb/mlhd/mRNA/2018-1-29_Myers/H_1105_BA9_mRNASeq__R1.fastq.gz', '/restricted/projectnb/mlhd/mRNA/2018-1-29_Myers/H_1105_BA9_mRNASeq__R2.fastq.gz'],
  ['C_0113_CAP_mRNASeq', '/restricted/projectnb/mlhd/mRNA/2018-1-29_Myers/C_0113_CAP_mRNASeq__R1.fastq.gz', '/restricted/projectnb/mlhd/mRNA/2018-1-29_Myers/C_0113_CAP_mRNASeq__R2.fastq.gz'],
  ['H_1104_CAP_mRNASeq', '/restricted/projectnb/mlhd/mRNA/2018-1-29_Myers/H_1104_CAP_mRNASeq__R1.fastq.gz', '/restricted/projectnb/mlhd/mRNA/2018-1-29_Myers/H_1104_CAP_mRNASeq__R2.fastq.gz'],
]

for dsid, r1, r2 in new_samples :

  fields = [ 'Dataset.dataset_id', 'mRNASeq.read1_path', 'mRNASeq.read2_path',
    'Dataset.datatype', 'Dataset.protocol', 'Dataset.batch', 'mRNASeq.readlen'
  ]
  vals = [
    dsid, r1, r2,
    'mRNASeq', 'TruSeq poly-A', 'mRNA_HD4', 101
  ]

  for f,v in zip(fields,vals) :
    df.loc[dsid,f] = v
 
  # remove and create the symlinks for these samples in the samples directory
  r1_link = '../samples/{}__R1.fastq.gz'.format(dsid)
  r2_link = '../samples/{}__R2.fastq.gz'.format(dsid)

  #try :
  #  os.remove(r1_link)
  #except :
  #  pass

  #try :
  #  os.remove(r2_link)
  #except :
  #  pass

  #os.symlink(r1, r1_link)
  #os.symlink(r2, r2_link)

df.to_csv('HD_mRNASeq_sample_info.csv',index=False)
