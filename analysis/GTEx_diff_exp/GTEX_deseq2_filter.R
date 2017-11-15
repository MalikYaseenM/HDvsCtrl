# DE Analysis of RNA-seq data
##load libraries
.libPaths( c( .libPaths(), '/projectnb/bubhub/conda_root/user_envs/crespodi/presymptomatic_hd_mrnaseq/lib/R/library') )
library(DESeq2)
## loading in the files
GTEX <- '../../samples/GTEx/results/GTEX_salmon_filter.csv'
meta_data <- '../../samples/GTEx/results/new_sample_info.csv'

##meta data
colData <- read.csv(meta_data, header=TRUE,stringsAsFactors=FALSE)
colData$SEX <- as.factor(colData$SEX)
colData$brain_region <- as.factor(colData$brain_region)

## GTEX data
GTEX_counts <- read.csv(GTEX, header=TRUE, row.names = 'gene_id')
GTEX_counts <- as.matrix(GTEX_counts)
storage.mode(GTEX_counts) <- 'integer'
GTEX_dds <- DESeqDataSetFromMatrix(countData = GTEX_counts,
                              colData = colData,
                              design = ~ AGE + SEX + brain_region)
GTEX_dds$brain_region<- relevel(GTEX_dds$brain_region, ref='CAU')
GTEX_dds <- DESeq(GTEX_dds)
GTEX_results <- results(GTEX_dds)
write.csv(GTEX_results, file = '../../samples/GTEx/results/GTEX_filter_deseq2_results.csv', quote = FALSE)
