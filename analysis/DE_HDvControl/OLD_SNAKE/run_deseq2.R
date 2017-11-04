library(docopt)
library(DESeq2)

'Usage: run_DEsep2.R <count_matrix>' -> doc
opts <- docopt(doc, commandArgs(trailingOnly = TRUE))

# read tsv file as csv and initialize counts as matrix
counts <- read.csv(opts$count_matrix, header=TRUE, row.names = 'gene_id')
counts <- as.matrix(counts)

# To prevent having DDS error of having missing integers use:
storage.mode(counts) <- 'integer'

# Taking out the sample names from the column names and saving it as a data frame
sample_names <- (colnames(counts))
colData <- as.data.frame(colnames(counts))
# Changing the row names into sample names
row.names(colData) <- colData$`colnames(counts)`

# Adding condition into the data frame to create colData
condition <- c()
for (i in 1:nrow(colData)) {
  sampleName <- row.names(colData)[i]
  if(grepl("^C", sampleName)) {
    condition <- c(condition,'control')
  } else {
    condition <- c(condition, 'HD')
  }
}

condition <- as.factor(condition)
colData$condition <- condition
# Deleting the extra column
colData$`colnames(counts)` <- NULL

# Setting the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ condition)
# Setting the reference to be control
dds$condition <- relevel(dds$condition, ref='control')
# Pre filtering rows 
dds <- dds[ rowSums(counts(dds)) > 1, ]
#dds <- dds[ rowMeans(counts(dds)) > 10, ]

# Deseq2 analysis
dds <- DESeq(dds)
res <- results(dds)
# Order results by smallest adjusted p value
resOrdered <- res[order(res$padj),]

# extracting normalized counts
dds <- estimateSizeFactors(dds)
# divides counts(dds) by sizeFactors(dds) by specifying normalized = TRUE
norm.counts <- counts(dds, normalized=TRUE)

name <- gsub('.tsv', '_normalized_counts.tsv',basename(opts$count_matrix))
colnames(norm.counts) <- sample_names
#write.table(norm.counts, file = file.path('samples',name), sep='\t', col.names = TRUE)
write.table(norm.counts, file = name, sep = ',', col.names = TRUE)


