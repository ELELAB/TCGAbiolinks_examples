####IMPORTANT! Please download the github version to use the most recent version of TCGAbiolinks
devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks")
###Use devtools::install_github(repo = "ELELAB/TCGAbiolinks")

library(TCGAbiolinks)
library(SummarizedExperiment)
library(recount)
###conversion of uuids to TCGA barcodes
library(TCGAutils)
library(limma)
library(biomaRt)


####Please download the github version to use the right version of TCGAbiolinks and TCGAutils

###Use devtools::install_github()


####Function to convert ENSG to Symbol
convert.ENSG.Symbol<-function(genes){
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
  #merge(DEG.ucs,G_list,by.x="gene",by.y="ensembl_gene_id")
  return(G_list)
}



########Query from Recount2 platform: Uterine Carcinoma#######
ucs.recount.gtex<-TCGAquery_recount2(project="GTEX", tissue="uterus")
ucs.recount.tcga<-TCGAquery_recount2(project="TCGA", tissue="uterus")

#### The steps below are needed to have the right correspondance beetween barcodes (TCGA) and UUID (recount)


query.ucs<- GDCquery(project = "TCGA-UCS",
                     data.category = "Transcriptome Profiling",
                     data.type = "Gene Expression Quantification",
                     workflow.type = "HTSeq - Counts")


samplesDown.ucs <- getResults(query.ucs,cols=c("cases"))


###tumor samples for uterine cancer
dataSmTP.ucs <- TCGAquery_SampleTypes(barcode = samplesDown.ucs,
                                      typesample = "TP")

###to check that there are no NT samples
dataSmNT.ucs <- TCGAquery_SampleTypes(barcode = samplesDown.ucs,
                                      typesample = "NT")


#####Preparing/scaling Recount2 data because it was sequenced using Rail-RNA

eset.gtex<-assays(scale_counts(ucs.recount.gtex$GTEX_uterus, round = TRUE))$counts
eset.tcga<-assays(scale_counts(ucs.recount.tcga$TCGA_uterus, round = TRUE))$counts

#### Check that the number of reads is less than or equal to 40 million
rse_scaled <- scale_counts(ucs.recount.gtex$GTEX_uterus, round = TRUE)
summary(colSums(assays(rse_scaled)$counts)) / 1e6

###replacing UUIDs with TCGA barcodes:
colnames(eset.tcga)<-colData(ucs.recount.tcga$TCGA_uterus)$gdc_cases.samples.portions.analytes.aliquots.submitter_id

###Removing version (number after ".")
rownames(eset.gtex) <- gsub("\\..*", "", rownames(eset.gtex))
rownames(eset.tcga) <- gsub("\\..*", "", rownames(eset.tcga))

####Segregate between primary tumors and normal samples
eset.tcga.cancer<-eset.tcga[,which(colData(ucs.recount.tcga$TCGA_uterus)$gdc_cases.samples.sample_type=="Primary Tumor")]
eset.tcga.nornal<-eset.tcga[,which(colData(ucs.recount.tcga$TCGA_uterus)$gdc_cases.samples.sample_type=="Solid Tissue Normal")]
####


##merging data by row names
dataPrep.ucs<-merge(as.data.frame(eset.gtex), as.data.frame(eset.tcga.cancer), by=0, all=TRUE)

rownames(dataPrep.ucs)<-dataPrep.ucs$Row.names
dataPrep.ucs$Row.names<-NULL


dataNorm.ucs <- TCGAanalyze_Normalization(tabDF = dataPrep.ucs,
                                          geneInfo = geneInfoHT,
                                          method = "gcContent")

dataFilt.ucs <- TCGAanalyze_Filtering(tabDF = dataNorm.ucs,
                                      method = "quantile", 
                                      qnt.cut =  0.25)

###set "metadata" argument to true when dealing with TCGA data
#in order to extract batch correction data
#processing 90 normal and 610 tumor samples on 17443 miRNA or genes
DEG.ucs <- TCGAanalyze_DEA( mat1 = dataFilt.ucs[,colnames(eset.gtex)],
                            mat2 = dataFilt.ucs[,colnames(eset.tcga.cancer)],
                            metadata =FALSE,
                            pipeline="limma",
                            voom = TRUE,
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")

View(DEG.ucs)

##These genes were removed because they cause duplicated rows
c("ENSG00000276085", "ENSG00000278233")
DEG.ucs<-DEG.ucs[-which(rownames(DEG.ucs) %in% c("ENSG00000276085", "ENSG00000278233")), ]

##converting ensenmbl gene ids to huugo gymbols using Biomart package
conversion.table<-convert.ENSG.Symbol(rownames(DEG.ucs))

conversion.inter.DEG<-intersect(conversion.table[-which(conversion.table$hgnc_symbol==""),]$ensembl_gene_id, rownames(DEG.ucs))
conversion.table2<-conversion.table[which(conversion.table$ensembl_gene_id %in% conversion.inter.DEG),]
rownames(conversion.table2)<-conversion.table2$ensembl_gene_id

conversion.table[-which(conversion.table$hgnc_symbol==""),]
DEG.ucs.hgnc<-DEG.ucs[conversion.inter.DEG,]

##after removing duplicate row names, we merge the hugo symboles with DEGs table
DEGs.ucs.hgnc<-merge(DEG.ucs.hgnc, conversion.table2, by=0)
rownames(DEGs.ucs.hgnc)<-DEGs.ucs.hgnc$Row.names
DEGs.ucs.hgnc$Row.names<-NULL

rownames(DEGs.ucs.hgnc)<-DEGs.ucs.hgnc$hgnc_symbol

write.csv(DEGs.ucs.hgnc, "DEGsUCS.csv")


###plot for genes with logfc > 5
TCGAVisualize_volcano(DEGs.ucs.hgnc$logFC, DEGs.ucs.hgnc$adj.P.Val,
                      filename = "UCS_+5.pdf", xlab = "logFC",
                      names = DEGs.ucs.hgnc$hgnc_symbol, show.names = "highlighted",
                      x.cut = 1, y.cut = 0.01, 
                      highlight =DEGs.ucs.hgnc$hgnc_symbol[which(DEGs.ucs.hgnc$logFC >= 5)][1:50],
                      highlight.color = "orange")

###plot for genes with logfc < -5
TCGAVisualize_volcano(DEGs.ucs.hgnc$logFC, DEGs.ucs.hgnc$adj.P.Val,
                      filename = "UCS_-5.pdf", xlab = "logFC",
                      names = DEGs.ucs.hgnc$hgnc_symbol, show.names = "highlighted",
                      x.cut = 1, y.cut = 0.01, 
                      highlight =DEGs.ucs.hgnc$hgnc_symbol[which(DEGs.ucs.hgnc$logFC <= -5)][1:50],
                      highlight.color = "orange")

