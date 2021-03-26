setwd("/Users/uqqzhao/Downloads")
count <- read.csv(file="gene_count_matrix.csv", head=T, sep=",")
se <- count

library("GenomicAlignments")
library("BiocParallel")


library("edgeR")
x <- read.csv(file="gene_count_matrix.csv", head=T, sep=",", row.names="gene_id")
group <- factor(c(2, 2, 2, 1, 1, 1))
y <- DGEList(counts=x,group=group)

keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

pdf("MDS.pdf")
col<-as.numeric(y$samples$group)
plotMDS(y,col=col)
dev.off()

design <- model.matrix(~group)
y <- estimateDisp(y,design)

fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
topTags(lrt)

summary(decideTests(lrt))

# Plot log-fold change against log-counts per million, with DE genes highlighted:
pdf("logFC_vs_logCPM.pdf")
plotMD(lrt)
abline(h=c(-1, 1), col="blue")
dev.off()

cpm <- cpm(y)
num_genes<-nrow(y$count)
diff<-topTags(lrt,n=num_genes)
# only take logFC,  
op1 <- diff$table[, c(1,4,5)]
merged <-  merge(op1, cpm, by=0, all=TRUE)
write.table(merged,"RCvsEXT.edgeR", row.names=FALSE, quote = FALSE, sep = "\t")

# save the sessionInfo
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

