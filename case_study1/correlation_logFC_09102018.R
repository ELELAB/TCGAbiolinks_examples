library(ggplot2)

DEGs_limma <- read.csv("DEGsBRCA_limma_091018.csv", row.names = 1)
DEGs_edgeR <- read.csv("DEGsBRCA_edgeR_091018.csv", row.names = 1)

genes_toTest <- rownames(DEGs_limma)[1:1000]
genes_common <- intersect(genes_toTest, rownames(DEGs_edgeR))

dataset <- subset(DEGs_limma, rownames(DEGs_limma) %in% genes_common, select = "Mycontrast.logFC")
dataset <- cbind(DEGs_edgeR$Mycontrast.logFC[match(rownames(dataset),rownames(DEGs_edgeR))], dataset)
colnames(dataset) <- c("logFC_edgeR", "logFC_limma")


pdf("scatterplot_logFC_limma_edgeR_top1000.pdf", width=5, height=3)
ggplot(dataset, aes(x=logFC_limma,y=logFC_edgeR))+geom_point()+
  xlab("logFC_limma")+ylab("logFC_edgeR")+
  theme(legend.title=element_blank(),axis.title=element_text(size=10),axis.text=element_text(size=10),
        legend.text=element_text(size=20))
dev.off()
cor(dataset$logFC_edgeR, dataset$logFC_limma)  #0.9936
