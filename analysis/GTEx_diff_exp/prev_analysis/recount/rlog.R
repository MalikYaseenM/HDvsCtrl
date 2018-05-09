library("DESeq2")
# filtered counts row means > 10
GTEX <-"../../samples/GTEx/results/GTEX_salmon_filter.csv"
cts <- as.matrix(read.csv(GTEX, header=TRUE, row.names = 'gene_id'))
storage.mode(cts)<- 'integer'

# perform rlog
rld = rlog(cts, blind=FALSE)

#write out to file
write.csv(rld, file="../../samples/GTEx/results/GTEX_salmon_rlog.csv")



