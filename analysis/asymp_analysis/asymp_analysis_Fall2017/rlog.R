library(DESeq2)
library(docopt)

#'Usage: Rscript rlog.R <count_matrix>' -> doc
#opts <- docopt(doc, commandArgs(trailingOnly = TRUE))

# Read salmon norm file
#norm <- read.csv(opts$count_matrix, header=T)
counts <- read.csv("../DE_HDvControl/all_salmon_norm.csv", header=T, row.names= "gene_id")
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
        condition <- c(condition,"control")
    } else {
        condition <- c(condition, "HD")
    }
}

condition <- as.factor(condition)
colData$condition <- condition
# Deleting the extra column
colData$`colnames(counts)` <- NULL

# Setting the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ condition)

# R log transformation
rld <- rlog(rld, blind=FALSE)
# Change matrix to data frame
rlog_output <- as.data.frame(assay(rld))
# Make new gene_id column
rlog_output$gene_id <- row.names(counts)
# Move gene_id to first
rlog_output <- rlog_output[,c(ncol(rlog_output),1:(ncol(rlog_output)-1))]

# Write csv output with row names as a new column
write.csv(rlog_output, "all_norm_rlog.csv", row.names=FALSE)
