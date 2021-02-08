library(xtail)

RNA_control1 <- 'jh46'
RNA_control2 <- 'jh63'
RNA_deg1 <- 'jh47'
RNA_deg2 <- 'jh64'
mrna <- read.csv("reads_FP_RNA.txt",row.names='gene', sep='\t')[ ,c(RNA_control1,RNA_control2,RNA_deg1,RNA_deg2)]
colnames(mrna) <- c("control1", "control2","deg1","deg2")

FP_control1 <- 'jh44'
FP_control2 <- 'jh61'
FP_deg1 <- 'jh45'
FP_deg2 <- 'jh62'
rpf <- read.csv("reads_FP_RNA.txt",row.names='gene', sep='\t')[ ,c(FP_control1,FP_control2,FP_deg1,FP_deg2)]
colnames(rpf) <- c("control1", "control2","deg1","deg2")

condition <- c("control", "control","deg","deg") 
test.results <- xtail(mrna,rpf,condition,bins=1000)

summary(test.results)

test.tab <- resultsTable(test.results)
write.csv(test.tab,"xtail_trmD_2way.csv",quote=F)