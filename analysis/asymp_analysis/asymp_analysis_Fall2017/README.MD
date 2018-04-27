Analysis between asymptomatic individuals (HDPos)

There are 5 asymptomatic individuals:
HD BA9: 'H_1105_BA9_mRNASeq', 'H_1104_BA9_mRNASeq', 'H_1106_BA9_mRNASeq'
HD CAP: 'H_1105_CAP_mRNASeq', 'H_1104_CAP_mRNASeq'

#### OLD ASYMP
Analysis steps (see Snakefile and python files for more details):
1. Filter 
2. Normalize
3. Firth
4. Analysis of normalized counts

Problem: 
1. Too many outliers -> did Spearman correlation to find out issue
2. Apparently, some of the samples were prepared differently and in different batches. 

Solution: Normalize using both batches.

Poly A: 'H_1105_BA9_mRNASeq', 'H_1104_BA9_mRNASeq'
Ribo Depleted: 'H_1106_BA9_mRNASeq', 'H_1105_CAP_mRNASeq', 'H_1104_CAP_mRNASeq'

#### NEW ASYMP
Analysis steps (see Snakefile and python files for more details):
1. Combine both batches:
Batch 1: 8 polyA-tail samples, 8 HD
Batch 2: 13 ribo-depleted samples, 7HD 6C
2. Filter 
3. Normalize
4. Pull out relevant samples
5. Firth
6. Analysis of normalized counts

Problem:
1. Still too many outliers
2. RNA processing batch effects can be seen from PCA of all samples 

#######################
