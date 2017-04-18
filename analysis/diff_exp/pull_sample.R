library(docopt)

'Usage: pull_sample.R <count_matrix>' -> doc
opts <- docopt(doc, commandArgs(trailingOnly = TRUE))

# read tsv file as csv and initialize counts as matrix
counts <- read.table(opts$count_matrix, header=TRUE, row.names = 'gene_id')
counts <- as.matrix(counts)

# Taking out the sample names from the column names and saving it as a data frame
colData <- as.data.frame(colnames(counts))

# Initialize ordered samples from colData
ordered <- colData

# Initialize original samples as ori
ori <- read.csv("../HD_mRNASeq_sample_info.csv", header=FALSE)

# Delete first row
ori <- ori[-1,]

# Set second row as column names
names(ori) <- as.matrix(ori[1, ])
ori <- ori[-1,]
ori[] <- lapply(ori,function(x) type.convert(as.character(x)))

# Make new matrix column15: datasetid, 2: subject type, 5: death age
design <- subset(ori, select = c(15, 2, 5))

# Set sample names as rownames
row.names(design) <- design$`Dataset.dataset_id`
row.names(ordered) <- ordered$`colnames(counts)`

# Finds similar row names and append subject type and death age to new data frame 
design <- design[row.names(ordered),]

# Write to new csv file
#write.table(design,"sample_info_design.csv",quote=FALSE,row.names=FALSE)

write.table(design,"sample_info_design.csv",col.names=TRUE,quote=FALSE,row.names=FALSE,sep=",")
