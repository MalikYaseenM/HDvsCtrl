# DE Analysis of RNA-seq data
##load libraries
.libPaths( c( .libPaths(), '/projectnb/bubhub/conda_root/user_envs/crespodi/presymptomatic_hd_mrnaseq/lib/R/library') )
library(DESeq2)
## loading in the files
GTEX <- '../../samples/GTEx/GTEx_salmon_filter.csv'
meta_data <- '../../samples/GTEx/GTEx_deseq_design.csv'

##meta data
colData <- read.csv(meta_data, header=TRUE,stringsAsFactors=FALSE, check.names=FALSE)
##colData$SEX <- as.factor(colData$SEX)
colData$brain_region <- as.factor(colData$brain_region)

## GTEX data
GTEX_counts <- read.csv(GTEX, header=TRUE, row.names='gene_id', check.names=FALSE)
## Get just GTEX names
GTEX_counts <- GTEX_counts[,grep("GTEX-", colnames(GTEX_counts))]
GTEX_counts <- as.matrix(GTEX_counts)
storage.mode(GTEX_counts) <- 'integer'
GTEX_dds <- DESeqDataSetFromMatrix(countData = GTEX_counts,
                              colData = colData,
                              design = ~ brain_region)
GTEX_dds$brain_region<- relevel(GTEX_dds$brain_region, ref='CAU')
GTEX_dds <- DESeq(GTEX_dds, minReplicatesForReplace = Inf)

## remove cooks cutoff outlier filtering
GTEX_results <- results(GTEX_dds, cooksCutoff = Inf)
write.csv(GTEX_results, file = '../../samples/GTEx/GTEX_filter_deseq2_results.csv', quote = FALSE)
