library(TCGAbiolinks)
library(plyr)
library(limma)
library(biomaRt)
library(SummarizedExperiment)

####IMPORTANT! Please download the github version to use the most recent version of TCGAbiolinks
###Use devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks")
###Use devtools::install_github(repo = "ELELAB/TCGAbiolinks")
################Querying, Downloading and batch correction#############################

query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")

samplesDown <- getResults(query,cols=c("cases"))

dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "TP")

dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "NT")

### for this example we will only use 100 samples for TP and 100  for NT ###
dataSmTP_short <- dataSmTP[1:100]
dataSmNT_short <- dataSmNT[1:100]

queryDown <- GDCquery(project = "TCGA-BRCA", 
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - Counts", 
                      barcode = c(dataSmTP_short, dataSmNT_short))

GDCdownload(query = queryDown)

dataPrep1 <- GDCprepare(query = queryDown, save = TRUE, save.filename = "brca_case1.rda")
#used GRCh38.p12 - From the 60483 genes we couldn't map 3767

#a step to remove sample outliers using pearson correlation
dataPrep <- TCGAanalyze_Preprocessing(object = dataPrep1, 
                                      cor.cut = 0.6,
                                      datatype = "HTSeq - Counts")

###############Tumor purity filtering###########

###vector containing all TCGA barcodes that hhave 60% tumor purity or more
Purity.BRCA<-TCGAtumor_purity(colnames(dataPrep1), 0, 0, 0, 0, 0.6)$pure_barcodes

################DEA with Molecular subtypes############

####Molecular subtypes:

####Filtering data so all samples have a pam50 subtype for BRCA:

##diff contains TCGA samples that have an available molecular subtype
###Also Applying Tumor purity filtering by making a set difference
#documentation for TCGA_MolecularSubtype is available

diff<-setdiff(Purity.BRCA, TCGA_MolecularSubtype(colnames(dataPrep[,dataSmTP_short]))$filtered)

dataFilt.brca.cancer<-dataPrep[,diff]
dataFilt.brca.normal<-dataPrep[,dataSmNT_short]
dataFilt.brca<-cbind(dataFilt.brca.normal, dataFilt.brca.cancer)


mol_subtypes<-c(rep("Normal", 100), TCGA_MolecularSubtype2(colnames(dataFilt.brca.cancer))$subtypes$subtype)
mol_subtypes<-make.names(mol_subtypes)

rownames(dataFilt.brca)<-rowData(dataPrep1)$external_gene_name


dataNorm.brca <- TCGAanalyze_Normalization(tabDF = dataFilt.brca,
                                           geneInfo = geneInfo,
                                           method = "gcContent")


dataFilt.brca.final <- TCGAanalyze_Filtering(tabDF = dataNorm.brca,
                                             method = "quantile", 
                                             qnt.cut =  0.25)



###Differential expression analysis after correcting for "TSS" factor.

DEG.brca.edgeR <- TCGAanalyze_DEA(MAT=dataFilt.brca.final,
                            pipeline="edgeR",
                            batch.factors = c("TSS"),
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            voom = FALSE,
                            method = "glmLRT", ClinicalDF = data.frame(),
                            Condtypes = mol_subtypes,
                            contrast.formula = "Mycontrast = (BRCA.LumA+BRCA.LumB)/2 -Normal")
#in this DEA we use Normal as a reference, thus genes with LogFC > 1 are up regulated in the subtypes
#with respect to the normal samples and viceversa for the LogFC < -1.

write.csv(DEG.brca.edgeR, "DEGsBRCA_edgeR_091018.csv", quote = FALSE)
#3404 genes identified

DEG.brca.edgeR.filt.TSS<-DEG.brca.edgeR$Mycontrast[which(abs(DEG.brca.edgeR$Mycontrast$logFC) >= 6), ]
#57 genes

TCGAVisualize_volcano(DEG.brca.edgeR$Mycontrast$logFC, DEG.brca.edgeR$Mycontrast$FDR,
                      filename = "LuminalABvsNormal_FC6.TSS.edgeR.pdf", xlab = "logFC",
                      names = rownames(DEG.brca.edgeR$Mycontrast), show.names = "highlighted",
                      x.cut = 1, y.cut = 0.01, 
                      highlight = rownames(DEG.brca.edgeR$Mycontrast)[which(abs(DEG.brca.edgeR$Mycontrast$logFC) >= 6)],
                      highlight.color = "orange")

