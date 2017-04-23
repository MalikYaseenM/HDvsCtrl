library(docopt)

'Usage: filter.R <count_matrix>' -> doc
opts <- docopt(doc, commandArgs(trailingOnly = TRUE))

# Read csv file as csv and initialize counts as matrix
counts <- read.csv(opts$count_matrix, header=TRUE, row.names = 'gene_id')
#counts <- as.matrix(counts)

# Initialize counts as a new variable
dfr <- counts

# Filter out rows with row means less than 10
dfr <- dfr[rowMeans(dfr) > 10,]

# Write to csv file
name <- gsub('quant_norm.csv', 'norm_filter.csv',basename(opts$count_matrix))
write.table(dfr, file = name, sep = ',', col.names = TRUE, quote=FALSE)
#write.table(head_all_salmon_norm_filter.csv, file = name, sep = ',', col.names = TRUE, quote=FALSE)
