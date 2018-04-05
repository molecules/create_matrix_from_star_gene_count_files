#!/bin/env Rscript
#SBATCH -J DESeq2 
#SBATCH -o out_DESeq2.o_%j
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=100G

library(edgeR)
library(DESeq2)

data = read.table("matrix.txt" header=T, row.names=1, com='')
rnaseqMatrix = round(data)
rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 1,]
conditions = data.frame(conditions=factor(c(rep("ForrM", 3), rep("ForrI", 3))))
rownames(conditions) = colnames(rnaseqMatrix)
ddsFullCountTable <- DESeqDataSetFromMatrix(
    countData = rnaseqMatrix,
    colData = conditions,
    design = ~ conditions)
dds = DESeq(ddsFullCountTable)
contrast=c("conditions","ForrM","ForrI")
res = results(dds, contrast)
baseMeanA <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$condition == "ForrM"])
baseMeanB <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$condition == "ForrI"])
res = cbind(baseMeanA, baseMeanB, as.data.frame(res))
res = cbind(sampleA="ForrM", sampleB="ForrI", as.data.frame(res))
res$padj[is.na(res$padj)]  <- 1
write.table(as.data.frame(res[order(res$pvalue),]), file='CountMatrix.txt.ForrM_vs_ForrI.DESeq2.DE_results', sep='  ', quote=FALSE)
write.table(rnaseqMatrix, file='CountMatrix.txt.ForrM_vs_ForrI.DESeq2.count_matrix', sep='  ', quote=FALSE)
source("/share/ircf/ircfapps/src/trinityrnaseq-Trinity-v2.3.2/Analysis/DifferentialExpression/R/rnaseq_plot_funcs.R")
pdf("MA_n_Volcano.pdf")
plot_MA_and_Volcano(log2(res$baseMean+1), res$log2FoldChange, res$padj)
dev.off()

