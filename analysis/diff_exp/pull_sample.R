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
write.table(sample_name,"sample_info.csv",col.names=FALSE,quote=FALSE,row.names=FALSE,sep=",")
