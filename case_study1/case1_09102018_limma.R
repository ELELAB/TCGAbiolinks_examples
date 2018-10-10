#setwd("~/test_ele2/")
library(TCGAbiolinks)
library(plyr)
library(limma)
library(biomaRt)
library(SummarizedExperiment)

####IMPORTANT! Please download the github version to use the most recent version of TCGAbiolinks
###Use devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks", force = TRUE)
###Use devtools::install_github(repo = "ELELAB/TCGAbiolinks")

####IMPORTANT2! If you want to reproduce exactly the data of this case study with the data downloaded at the 
#### time of the analysis you would need to comment the GDCquery and
####GDCprepare steps and starting loading the .rda file already assembled by us using
####dataPrep1 <- get(load("brca_case1.rda"))


################Querying, Downloading and batch correction#############################

query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")

samplesDown <- getResults(query,cols=c("cases"))
#1222 samples
dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "TP")
#1102
dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "NT")
#113

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

#step with library size and gcContent normalization using EDASeq
dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfoHT,
                                      method = "gcContent")
#quantile filtering to remove genes with low count
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                    qnt.cut =  0.25)
#voom transformation of the data (log)
v.dataFilt<-voom(dataFilt)

#in order to have 2 batches with multiple samples and avoid batches with one sample
#the user need to call first the get_IDs function on the top of this script if this has not been done already
c1<-which(get_IDs(dataPrep)$tss=="E9")
c2<-which(get_IDs(dataPrep)$tss=="E2")

#taking log transformed data for exploration of batch effects

TCGAbatch_Correction(tabDF = v.dataFilt$E[,c(c1,c2)], batch.factor="TSS", adjustment=NULL)
#if you get error in plot.new():figure margins too large - please use:
#par(mar=c(1,1,1,1))
#and then rerun the command above

#### After exploration of the batches we can continue to work on the original gene matrix for DEA ###
###############Tumor purity filtering###########

###vector containing all TCGA barcodes that hhave 60% tumor purity or more
Purity.BRCA<-TCGAtumor_purity(colnames(dataPrep1), 0, 0, 0, 0, 0.6)$pure_barcodes

################DEA with Molecular subtypes############

####Molecular subtypes:

####Filtering data so all samples have a pam50 subtype for BRCA:

##diff contains TCGA samples that have an available molecular subtype
###Also Applying Tumor purity filtering by making a set difference
#documentation for TCGA_MolecularSubtype is available

##########Important########
####in this example no change happens because all the TCGA barcodes
####have available subtype info

###setdiff() is here used to remove samples with no subtype info from the TCGA BRCA samples that satisfy
###tumor purity criteria

diff<-setdiff(Purity.BRCA, TCGA_MolecularSubtype(colnames(dataPrep[,dataSmTP_short]))$filtered)

dataFilt.brca.cancer<-dataPrep[,diff]
dataFilt.brca.normal<-dataPrep[,dataSmNT_short]
dataFilt.brca<-cbind(dataFilt.brca.normal, dataFilt.brca.cancer)


mol_subtypes<-c(rep("Normal", 100), TCGA_MolecularSubtype(colnames(dataFilt.brca.cancer))$subtypes$subtype)
#All barcodes have available molecular subtype info
mol_subtypes<-make.names(mol_subtypes)

#to convert ENSEMBL id to gene name
rownames(dataFilt.brca)<-rowData(dataPrep1)$external_gene_name

#we redo here the normalization and filtering steps 
dataNorm.brca <- TCGAanalyze_Normalization(tabDF = dataFilt.brca,
                                           geneInfo = geneInfo,
                                           method = "gcContent")


dataFilt.brca.final <- TCGAanalyze_Filtering(tabDF = dataNorm.brca,
                                       method = "quantile", 
                                       qnt.cut =  0.25)


###Differential expression analysis after correcting for "TSS" factor.

DEG.brca.TSS <- TCGAanalyze_DEA(MAT=dataFilt.brca.final,
                            pipeline="limma",
                            batch.factors = c("TSS"),
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            voom = TRUE,
                            Condtypes = mol_subtypes,
                            contrast.formula = "Mycontrast = (BRCA.LumA+BRCA.LumB)/2 -Normal")
#in this DEA we use Normal as a reference, thus genes with LogFC > 1 are up regulated in the subtypes
#with respect to the normal samples and viceversa for the LogFC < -1.

write.csv(DEG.brca.TSS, "DEGsBRCA_limma_091018.csv", quote = FALSE)
#3241 genes identified

#to check how many with logFC > 5 (for plotting purposes)
DEG.brca.filt.TSS<-DEG.brca.TSS$Mycontrast[which(abs(DEG.brca.TSS$Mycontrast$logFC) >= 6), ]
#47 genes identified

TCGAVisualize_volcano(DEG.brca.TSS$Mycontrast$logFC, DEG.brca.TSS$Mycontrast$adj.P.Val,
                      filename = "LuminalABvsNormal_FC6.TSS.pdf", xlab = "logFC",
                      names = rownames(DEG.brca.TSS$Mycontrast), show.names = "highlighted",
                      x.cut = 1, y.cut = 0.01, 
                      highlight = rownames(DEG.brca.TSS$Mycontrast)[which(abs(DEG.brca.TSS$Mycontrast$logFC) >= 6)],
                      highlight.color = "orange")
