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

dataPrep1 <- GDCprepare(query = queryDown)

dataPrep <- TCGAanalyze_Preprocessing(object = dataPrep1, 
                                      cor.cut = 0.6,
                                      datatype = "HTSeq - Counts")

dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfoHT,
                                      method = "gcContent")


dataNorm <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                    qnt.cut =  0.25)

v.dataNorm<-voom(dataNorm)

#in order to have 2 batches with multiple samples and avoid batches with one sample
#the user need to call first the get_IDs function on the top of this script if this has not been done already
c1<-which(get_IDs(dataPrep)$tss=="E9")
c2<-which(get_IDs(dataPrep)$tss=="E2")

#taking log transformed data for exploration of batch effects

TCGAbatch_Correction(tabDF = v.dataNorm$E[,c(c1,c2)], batch.factor="TSS", adjustment=NULL)


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


mol_subtypes<-c(rep("Normal", 100), TCGA_MolecularSubtype(colnames(dataFilt.brca.cancer))$subtypes$subtype)
mol_subtypes<-make.names(mol_subtypes)

rownames(dataFilt.brca)<-rowData(dataPrep1)$external_gene_name


dataNorm.brca <- TCGAanalyze_Normalization(tabDF = dataFilt.brca,
                                           geneInfo = geneInfo,
                                           method = "gcContent")


dataFilt.brca.final <- TCGAanalyze_Filtering(tabDF = dataNorm.brca,
                                       method = "quantile", 
                                       qnt.cut =  0.25)

###Differential expression analysis after correcting for "Plate" factor.

DEG.brca <- TCGAanalyze_DEA(MAT=dataFilt.brca.final,
                            pipeline="limma",
                            batch.factors = c("Plate"),
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            voom = TRUE,
                            method = "glmLRT", ClinicalDF = data.frame(),
                            Condtypes = mol_subtypes,
                            contrast.formula = "Mycontrast = Normal - (BRCA.Luminal.A+BRCA.Luminal.B)/2")

write.csv(DEG.brca, "DEGsBRCA.csv")

DEG.brca.filt<-DEG.brca$Mycontrast[which(abs(DEG.brca$Mycontrast$logFC) >= 5), ]

TCGAVisualize_volcano(DEG.brca$Mycontrast$logFC, DEG.brca$Mycontrast$adj.P.Val,
                      filename = "LuminalABvsNormnal_5.pdf", xlab = "logFC",
                      names = rownames(DEG.brca$Mycontrast), show.names = "highlighted",
                      x.cut = 1, y.cut = 0.01, 
                      highlight = rownames(DEG.brca$Mycontrast)[which(abs(DEG.brca$Mycontrast$logFC) >= 5)],
                      highlight.color = "orange")


