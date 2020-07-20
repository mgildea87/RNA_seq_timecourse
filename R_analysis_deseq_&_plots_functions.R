setwd('Z:/mag456/Srb2/RNA-seq/Class Time course Fall 2017/RNA-seq analysis/')

#read in data
raw_cere_counts = read.table('raw_counts_cerevisiae.txt')
raw_pombe_counts = read.table('raw_counts_pombe.txt')
raw_total_counts = rbind(raw_cere_counts, raw_pombe_counts)


#draw MA plots

plotma = function(a ,b){
  library(ggplot2)
  library(zoo)
  library(lattice)
  M = vector(length = nrow(raw_total_counts))
  fc = vector(length = nrow(raw_total_counts))
  for (i in 1:nrow(raw_total_counts)){
    M[i] = log2(mean(c(raw_total_counts[i,a], raw_total_counts[i,b])))
    fc[i] = log2(raw_total_counts[i,a]/raw_total_counts[i,b])
  }
  data = data.frame(cbind(M, fc))
  mv = rollapply(data[,2],10, FUN=function(x) mean(x, na.rm = TRUE), partial = TRUE)
  ggplot(data = data, aes(x = data[,1], y = data[,2])) +
    geom_point(color = "blue", size = 0.75)+
    geom_point(data = data[6573:11690,], aes(x = data[6573:11690,1], y = data[6573:11690,2]), color = "red",
               size = 0.75)+
    geom_smooth(data = data[6573:11690,],aes(data[6573:11690,1], data[6573:11690,2]),method = "loess", color = "red3") +
    geom_hline(yintercept = 0, color = "grey") +
    geom_smooth(data = data[1:6573,],aes(data[1:6573,1], data[1:6573,2]),method = "loess", color = "blue3")+
    theme(panel.background = element_rect(fill = "grey90"))+
    labs(x = "Log2(mean(rep_A, rep_B))", y = "Log2 Fold Change A/B")
    
  
}  

#Diff exp following RUVseq documentation with DESEQ2
library(RUVSeq)
library(DESeq2)
library(Rcpp)
library(vsn)
library(RColorBrewer)
library(gplots)
data("pasillaGenes")
countData <- counts(pasillaGenes)
colData <- pData(pasillaGenes)[,c("condition","type")]

#prepare sample information data frame
Sample_data_mat <- matrix(nrow = ncol(raw_cere_counts), ncol = 1)
row.names(Sample_data_mat) <- colnames(raw_cere_counts)
  for (i in 1:nrow(Sample_data_mat)){
    x = row.names(Sample_data_mat)[i]
    Sample_data_mat[i,1] = substr(x, 1, nchar(x)-5)
  }

Sample_data = as.data.frame(Sample_data_mat)
colnames(Sample_data) <- "condition"

#correct raw pombe counts
correct_index = c(12, 14, 15, 16, 17, 18, 19, 20, 21, 23, 25, 27 , 28, 30, 32, 34, 36, 37)
for (i in correct_index){
  for (j in 1:nrow(raw_pombe_counts)){
    raw_pombe_counts[j,i] = raw_pombe_counts[j,i]/2
  }
}
#DESeq2
deseq = function(a,b){
dds <- DESeqDataSetFromMatrix(countData = raw_cere_counts, colData = Sample_data, design = ~ condition)
dds <- dds[, dds$condition %in% c(a,b) ]
dds$condition <- droplevels(dds$condition)
first = Sample_data$condition[min(which(a == Sample_data$condition)[1], which(b == Sample_data$condition)[1])]
second = Sample_data$condition[max(which(a == Sample_data$condition)[1], which(b == Sample_data$condition)[1])]
sfs <- estimateSizeFactorsForMatrix(raw_pombe_counts)
sfs <- sfs[c(which(Sample_data$condition == first),which(Sample_data$condition == second))]
sizeFactors(dds) = sfs
sizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
res <- results(dds)
res <- res[order(res$padj),]
write.table(res, paste("DESeq2_results/",a,"_vs_",b,".txt", sep = ""))

rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:100]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

par(mfrow=c(1,3))
notAllZero <- (rowSums(counts(dds))>0)
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
pdf(paste("DESeq2_results/",a,"_vs_",b, ".pdf", sep = ""), height = 5) 
plotMA(res)
meanSdPlot(log2(counts(dds,normalized=TRUE)[notAllZero,] + 1),ylim = c(0,2.5))
meanSdPlot(assay(rld[notAllZero,]),ylim = c(0,2.5))
meanSdPlot(assay(vsd[notAllZero,]),ylim = c(0,2.5))
heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10,6))
heatmap.2(assay(rld)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))
heatmap.2(assay(vsd)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))
plotDispEsts(dds)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
plotPCA(rld)
dev.off()
}
for (i in Sample_data$condition){
  for (j in Sample_data$condition){
    if (i == j){
      next
    }
    else{
    deseq(i, j)
    }  
  }
}


  
