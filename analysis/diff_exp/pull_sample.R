library(docopt)

'Usage: pull_sample.R <count_matrix>' -> doc
opts <- docopt(doc, commandArgs(trailingOnly = TRUE))

# read tsv file as csv and initialize counts as matrix
counts <- read.table(opts$count_matrix, header=TRUE, row.names = 'gene_id')
counts <- as.matrix(counts)

# Taking out the sample names from the column names and saving it as a data frame
sample_names <- (colnames(counts))
colData <- as.data.frame(colnames(counts))

# Make new colData as sample names
# Some colData had __salmon__counts ending (on sample names with REPEAT)
sample_name <- sub("__salmon__counts","",colData[,1])

# Sample names didn't have mRNASeq, added mRNASeq to all
sample_name <- sub("CAP$","CAP_mRNASeq",sample_name)
sample_name <- sub("BA9$","BA9_mRNASeq",sample_name)

# Set sample name as a data frame
sample_name <- as.data.frame(sample_name)

# Write to csv file
write.table(sample_name,"sample_info_norm_order.csv",col.names=FALSE,quote=FALSE,row.names=FALSE,sep=",")

# Initialize ordered samples as ordered
ordered <- read.csv("sample_info_norm_order.csv", header=FALSE)
ordered <- as.data.frame(ordered)

# Initialize original samples from python script output
original <- read.csv("sample_info_py.csv", header=FALSE)
original <- as.data.frame(original)
row.names(original) <- original$V1
row.names(ordered) <- ordered$V1

# Reorder original (append into ordered) and rename as original2
original2 <- original[row.names(ordered),]

# Pull subject type and age of death from original2
ordered_design <- original2[,(2:3)]
# Merge the 2 columns and add ~
ordered_design2 <- as.data.frame(paste(original2$V2, original2$V3, sep=" ~ "))

# Write to csv file
write.table(ordered_design2,"sample_info_design.csv",col.names=FALSE,quote=FALSE,row.names=FALSE,sep=",")
