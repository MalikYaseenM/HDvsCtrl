#library(docopt)
library(DESeq2)

#'Usage: run_DEseq2.R <count_matrix> <design_info>' -> doc
#opts <- docopt(doc, commandArgs(trailingOnly = TRUE))

# read csv file as csv and initialize counts as matrix
#counts <- read.csv(opts$count_matrix, header=TRUE, row.names = 'gene_id')
#samp <- read.csv(opts$design_info, header=TRUE, sep = ",")

count_matrix <- '../../samples/Analysis_Results/CAP_raw_filter.csv'
counts <- read.csv(count_matrix, header=TRUE, row.names = 'gene_id')
counts <- as.matrix(counts)

des_info <- '../HD_mRNASeq_sample_info.csv'
samp <- read.csv(des_info, header=TRUE, sep = ",")

# To prevent having DDS error of having missing integers use:
storage.mode(counts) <- 'integer'

# Taking out the sample names from the column names and saving it as a data frame)
colData <- as.data.frame(colnames(counts))
# Changing the row names into sample names
row.names(colData) <- colData$`colnames(counts)`

# Adding condition into the data frame to create colData
des <- subset(samp, select = c(15, 2, 5))
row.names(des) <- des$`Dataset.dataset_id`
colData <- des[match(colData$`colnames(counts)`, des$Dataset.dataset_id), ]
colData$Subject.subject_type <- as.factor(colData$Subject.subject_type)

# Deleting the extra column
colData$`colnames(counts)` <- NULL

# Setting the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ Subject.death + Subject.subject_type)
# Setting the reference to be control
dds$Subject.subject_type <- relevel(dds$Subject.subject_type, ref='Control')

# Deseq2 analysis
dds <- DESeq(dds, minReplicatesForReplace = Inf)
dds.results <- results(dds, cooksCutoff = Inf, independentFiltering=FALSE)
# Order results by smallest adjusted p value
#resOrdered <- res[order(res$padj),]

name <- gsub('.csv', '_deseq2.csv', count_matrix)
#write.table(norm.counts, file = file.path('samples',name), sep='\t', col.names = TRUE)
#write.table(norm.counts, file = name, sep = ',', col.names = TRUE)
write.table(dds.results, file = name, sep = ',')
