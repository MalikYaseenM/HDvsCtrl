library('recount')
getwd()
load('rse_gene_brain.Rdata')
##load library "library('recount')
##>rse_gene tells you about the rdata object
##>colnames(rse_gene) gives you the column names
##sampid contains the sample ids in the format that we know
##assay(1) gives you information about the data
##assay(rse_gene, "counts")[1:4, 1:4] first 4 rows and first 4 samples
##each row has an associated grange with it



## class: RangedSummarizedExperiment 
## dim: 58037 1409  58k rolls and 1409 samples
## metadata(0):
## assays(1): counts
## rownames(58037): ENSG00000000003.14 ENSG00000000005.5 ...
##   ENSG00000283698.1 ENSG00000283699.1
## rowData names(3): gene_id bp_length symbol
## colnames(1409): SRR2166176 SRR2167642 ... SRR1073143 SRR627421
## colData names(82): project sample ... title characteristics




                                        #url <- download_study('SRP025982', type='counts-gene')
## 




                                        #load(file.path('SRP025982', 'counts-gene'))

##sample info in sampid column. rse_gene$sampid
## assay(rse_gene, "counts")[1:4, 1:4]


assay(rse_gene, "counts")[1:4, 1:4]
total_experiment = assay(rse_gene, "counts")
write.csv(total_experiment, file="RSE_Counts.csv", row.names=TRUE)
