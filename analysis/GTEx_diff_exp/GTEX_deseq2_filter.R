# DE Analysis of RNA-seq data
##load libraries
.libPaths( c( .libPaths(), '/projectnb/bubhub/conda_root/user_envs/crespodi/presymptomatic_hd_mrnaseq/lib/R/library') )
library(DESeq2)

GTEX <- '../../samples/GTEx/GTEx_salmon_filter.csv'
meta_data <- '../../samples/GTEx/GTEx_deseq_design.csv'

GTEX_counts <- read.csv(GTEX, header=TRUE, row.names='gene_id', check.names=FALSE)

## Get just GTEX names
GTEX_counts <- GTEX_counts[,grep("GTEX-", colnames(GTEX_counts))]
GTEX_counts <- as.matrix(GTEX_counts)
storage.mode(GTEX_counts) <- 'integer'

meta_df <- read.csv(meta_data, header=TRUE,stringsAsFactors=FALSE, check.names=FALSE)
# Taking out the sample names from the column names and saving it as a data frame
colData <- as.data.frame(colnames(GTEX_counts))
# Changing the row names into sample names
row.names(colData) <- colData$`colnames(GTEX_counts)`

row.names(meta_df) <- meta_df$`identifier`
colData <- meta_df[match(colData$`colnames(GTEX_counts)`, meta_df$identifier), ]

colData$brain_region <- as.factor(colData$brain_region)

# Tests to check if the row and column names match in same order
all(rownames(colData) %in% colnames(GTEX_counts))
all(rownames(colData) == colnames(GTEX_counts))


GTEX_dds <- DESeqDataSetFromMatrix(countData = GTEX_counts,
                              colData = colData,
                              design = ~ brain_region)

GTEX_dds$brain_region<- relevel(GTEX_dds$brain_region, ref='BA9')

GTEX_dds <- DESeq(GTEX_dds, minReplicatesForReplace = Inf)

## remove cooks cutoff outlier filtering
GTEX_results <- results(GTEX_dds, cooksCutoff = Inf)
write.csv(GTEX_results, file = '../../samples/GTEx/GTEX_filter_deseq2_results.csv', quote = FALSE)




