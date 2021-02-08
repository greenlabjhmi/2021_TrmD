library(DESeq2)

all_counts = read.csv("reads_rnaseq.txt", row.names='gene', sep='\t')
all_counts = round(all_counts)


#import experimental design
coldata <- read.csv('deseq_design.txt', row.names='name', sep='\t')
coldata$genotype <- factor(coldata$genotype)

#test that all rows and columns in right order
all(rownames(coldata) == colnames(all_counts))

#define DESeq dataset 
dds <- DESeqDataSetFromMatrix(countData=all_counts, colData=coldata, design=~genotype)
dds$genotype <- relevel(dds$genotype, ref = "WT")

#run differential expression test
dds <- DESeq(dds)
#names of the resulting comparisons
resultsNames(dds)

res <-results(dds, contrast=c("genotype","mutant","WT"))
write.table(res,"mutant_vs_WT.tsv", quote=F, sep='\t', col.names=NA)
plotMA(res, ylim=c(-5,5))