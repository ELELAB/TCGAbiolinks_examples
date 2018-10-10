library(UpSetR)

###### plot the overlap between DEGs detected in limma and edgeR #########

DEA_limma <- read.csv("DEGsBRCA_limma_091018.csv", row.names = 1)
DEA_edgeR <- read.csv("DEGsBRCA_edgeR_091018.csv", row.names = 1)

up_limma <- rownames(DEA_limma)[which(DEA_limma$Mycontrast.logFC > 1 & DEA_limma$Mycontrast.adj.P.Val <0.01)]
down_limma <- rownames(DEA_limma)[which(DEA_limma$Mycontrast.logFC < -1 & DEA_limma$Mycontrast.adj.P.Val <0.01)]

up_edgeR <- rownames(DEA_edgeR)[which(DEA_edgeR$Mycontrast.logFC > 1 & DEA_edgeR$Mycontrast.FDR < 0.01)]
down_edgeR <- rownames(DEA_edgeR)[which(DEA_edgeR$Mycontrast.logFC < -1 & DEA_edgeR$Mycontrast.FDR < 0.01)]


listInput <- list(up_limma,up_edgeR,down_limma,down_edgeR)
names(listInput) <- c("up_limma","up_edgeR","down_limma", "down_edgeR")

pdf("overlap_limma-edgeR.pdf",width=10, height=7, onefile = FALSE)
upset(fromList(listInput), order.by = "freq",text.scale = c(2,3,2,1.8,3,3),
      mainbar.y.label = "DEGs Intersection", sets.x.label = "DEGs per method",
      line.size = 1)
dev.off()



