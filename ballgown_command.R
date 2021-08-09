if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ballgown")

library(ballgown)

library(dplyr) # for the arrange function

bg = ballgown(samples=c("RC1.ctab", "RC2.ctab", "RC3.ctab", "EXT1.ctab", "EXT2.ctab", "EXT3.ctab"), meas='all')
pData(bg) = data.frame(id=sampleNames(bg), group=c(0,0,0,1,1,1))

stat_results = stattest(bg, feature='transcript', meas='FPKM', getFC=TRUE, covariate='group')

results_transcripts = data.frame(transcriptIDs=ballgown::transcriptNames(bg), geneIDs=ballgown::geneIDs(bg), geneNames=ballgown::geneNames(bg), stat_results[,3:5], texpr(bg, 'FPKM'))

results_transcripts = arrange(results_transcripts,pval)
#head(results_transcripts)

write.table(results_transcripts,"Nucleus.lncRNACapture.ballgown.xls", row.names=FALSE, quote = FALSE, sep = "\t")
