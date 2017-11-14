# DE Analysis of RNA-seq data

##load libraries
## design=~ SEX + AGE 
.libPaths( c( .libPaths(), '/projectnb/bubhub/conda_root/user_envs/crespodi/presymptomatic_hd_mrnaseq/lib/R/library') )
library(DESeq2)
## loading in the files
RSE_file <- '../../samples/GTEx/results/filtered_RSE.csv'
meta_data <- '../../samples/GTEx/results/new_sample_info.csv'
##meta data
colData <- read.csv(meta_data, header=TRUE,stringsAsFactors=FALSE)
colData$SEX <- as.factor(colData$SEX)
colData$brain_region <- as.factor(colData$brain_region)

##RSE data
counts <- read.csv(RSE_file, header=TRUE, row.names = 'gene_id')
counts <- as.matrix(counts)
storage.mode(counts) <- 'integer'
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ AGE + SEX)
dds$brain_region<- relevel(dds$brain_region, ref='CAU')
dds <- DESeq(dds)
res <- results(dds)
write.csv(res, file = '../../samples/GTEx/results/NEW_RSE_filter_deseq2_results.csv', quote = FALSE)
