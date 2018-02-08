# DE Analysis of RNA-seq data

##load libraries
## design=~ SEX + AGE 
.libPaths( c( .libPaths(), '/projectnb/bubhub/conda_root/user_envs/crespodi/presymptomatic_hd_mrnaseq/lib/R/library') )
library(DESeq2)
## loading in the files
combat_file <- '../../samples/GTEx/results/GTEx_combat_newcol.csv'
meta_data <- '../../samples/GTEx/results/sample_info_design.csv'

##meta data
colData <- read.csv(meta_data, header=TRUE, stringsAsFactors=FALSE)
colData$SEX <- as.factor(colData$SEX)
colData$AGE <- as.factor(colData$AGE)
colData$brain_region <- as.factor(colData$brain_region)

##RSE data
counts <- read.csv(combat_file, header=TRUE, row.names = 'gene_id')
counts <- as.matrix(counts)
storage.mode(counts) <- 'integer'
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ brain_region)
dds$brain_region<- relevel(dds$brain_region, ref='CAU')
dds <- DESeq(dds, minReplicatesForReplace = Inf)
res <- results(dds, cooksCutoff = Inf)
write.csv(res, file = '../../samples/GTEx/results/GTEX_covar_deseq2_results.csv', quote = FALSE)
