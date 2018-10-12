##setwd("~/ele_test3/")
library(TCGAbiolinks)
library(plyr)
library(SummarizedExperiment)
library(limma)

#retrieve TCGA LUAD samples
query <- GDCquery(project = "TCGA-LUAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")

samplesDown <- getResults(query,cols=c("cases"))

dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "TP")

dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "NT")

queryDown <- GDCquery(project = "TCGA-LUAD", 
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - Counts", 
                      barcode = c(dataSmTP))

GDCdownload(query = queryDown)

dataPrep1 <- GDCprepare(query = queryDown)

#pre-processing steps for the TCGA data (not needed for GEO in this
#specific case of study since they have been already carried out by the authors of the GEO dataset
dataPrep <- TCGAanalyze_Preprocessing(object = dataPrep1, 
                                      cor.cut = 0.6,
                                      datatype = "HTSeq - Counts")

dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfoHT,
                                      method = "gcContent")

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)
#if you want to reproduce exactly this case of study we reccomand to start from here loading the dataFilt_LUAD.rda
#object that we provide in the example folder

save(dataFilt, file = "dataFilt_LUAD_tp.rda")

# retrieve only tumor samples
sampleTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), typesample = "TP")
dataFiltTP <- dataFilt[,sampleTP[1:100]]

# retrieve year of diagnosis
clin_LUAD <- GDCquery_clinic("TCGA-LUAD", type = "clinical", save.csv = TRUE)
clin_LUAD$age_at_diag_year <- floor(clin_LUAD$age_at_diagnosis/365)
clin_LUAD$diag_year <- clin_LUAD$age_at_diag_year + clin_LUAD$year_of_birth
# select only samples taken between 2010 and 2013
clin_LUAD <- subset(clin_LUAD, clin_LUAD$diag_year>=2010 & clin_LUAD$diag_year<=2012)
clin_LUAD <- clin_LUAD[which(clin_LUAD$submitter_id %in% substr(colnames(dataFiltTP),1,12)),]
dataFiltTP <- dataFiltTP[,which(substr(colnames(dataFiltTP),1,12) %in% clin_LUAD$submitter_id)]
dim(dataFiltTP)

#convert ensembl ID in gene name
rownames(dataFiltTP) <- rowData(dataPrep1)$external_gene_name[match(rownames(dataFiltTP),rowData(dataPrep1)$ensembl_gene_id)]

# log2-transgorm TCGA data to make them comparable to the GEO data
dataFiltTP <- log2(dataFiltTP[,]+1)

####################### GEO datset: GSE60052 ##############################

clin_GEO <- read.csv("GSE60052_clinical.csv", header = T, row.names = 1, check.names = F)

# retrieve only the year from whole date
date <- strsplit(as.character(clin_GEO$`Dx Date`),"/")
date <- ldply(date, rbind)
clin_GEO$year <- as.numeric(paste(date$`3`))
clin_GEO <- subset(clin_GEO, select = "year")

# select only samples taken between 2010 and 2013

clin_GEO <- subset(clin_GEO, clin_GEO$year>=2010 & clin_GEO<=2012)
dataset <- read.csv("GSE60052_79tumor.7normal.normalized.log2.data.Rda.tsv", sep = "\t", header = T, row.names = 1, check.names = F)

# remove normal samples
dataset <- dataset[,-grep("normal", colnames(dataset))]
# convert the sample names in order to be matched with that ones in clinical data
colnames(dataset) <- gsub("[[:alpha:]]", "", colnames(dataset))

#remove duplicates
dataset <- dataset[, -which(duplicated(colnames(dataset)) | duplicated(colnames(dataset), fromLast = T))]
dataset <- dataset[,colnames(dataset) %in% rownames(clin_GEO)]

# retrieve the genes in common between GEO and TCGA-LUAD datasets
dataset <- dataset[rownames(dataset) %in% intersect(rownames(dataset),rownames(dataFiltTP)),]
dataFiltTP <- dataFiltTP[rownames(dataFiltTP) %in% intersect(rownames(dataFiltTP),rownames(dataset)),]

# merge the two counts matrices
countsTable <- cbind(dataFiltTP,dataset[match(rownames(dataFiltTP), rownames(dataset)),])

############ pca ###################

pca <- prcomp(t(countsTable), center = TRUE, scale=TRUE)
pca <- as.data.frame(pca$x)

my.group <- rep(NA,ncol(countsTable))
my.group[which(colnames(countsTable) %in% TCGAquery_SampleTypes(barcode = colnames(countsTable), typesample = "TP"))] <- "TCGA.tumor"
my.group[which(is.na(my.group)==TRUE)] <- "GEO.tumor"
pca$group <- as.factor(my.group)

library (ggplot2)

pdf("pca_GEO_TCGA_log.pdf")
p<-ggplot(pca,aes(x=PC1,y=PC2,color=group ))
p<-p+geom_point()
print(p)
dev.off()

###################################################################

#create dataframe with batch information
AnnotationCounts <- matrix(0,ncol(countsTable),2)
colnames(AnnotationCounts) <- c("Samples","Batch")
rownames(AnnotationCounts) <- colnames(countsTable)
AnnotationCounts <- as.data.frame(AnnotationCounts)
AnnotationCounts$Samples <- colnames(countsTable)
AnnotationCounts[colnames(dataFiltTP),"Batch"] <- clin_LUAD$diag_year[match(substr(colnames(dataFiltTP),1,12), clin_LUAD$submitter_id)]
AnnotationCounts[colnames(dataset),"Batch"] <- clin_GEO$year[match(colnames(dataset),rownames(clin_GEO))]


countsCorrected <- TCGAbatch_Correction(tabDF = countsTable,
                                        UnpublishedData = TRUE, 
                                        AnnotationDF = AnnotationCounts)

write.csv(countsCorrected, "countsCorrected_GEO-TCGA.csv")

sampleTP <- TCGAquery_SampleTypes(barcode = colnames(countsCorrected), typesample = "TP")
tumorGEO <- colnames(countsCorrected)[which(!colnames(countsCorrected) %in% sampleTP==TRUE)]

DEGs <- TCGAanalyze_DEA(mat1 = countsCorrected[,tumorGEO],
                        mat2 = countsCorrected[,sampleTP],
                        Cond1type = "tumorGEO",
                        Cond2type = "tumorTCGA",
                        pipeline = "limma",
                        voom = FALSE,
                        fdr.cut = 0.005,
                        logFC.cut = 2)

TCGAVisualize_volcano(DEGs$logFC, DEGs$adj.P.Val,
                      filename = "TCGA-GEO.pdf", xlab = "logFC",
                      x.cut = 2, y.cut = 0.005)
