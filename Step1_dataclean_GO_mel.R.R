#DEGs claenup for enricgGO, gseGO and gseKEGG
#ONLY for melanogater use 

library(ggplot2)
library(DOSE)
library(pathview)
library("AnnotationDbi")
library("org.Dm.eg.db")
library("ggnewscale")
library(GOSemSim)
library(enrichplot)
library(clusterProfiler)

#read the file with has DEGs at padj < 0.05 FoldChange >=2
#Change the name based on the info in the .csv file 
#one of the col should contain the ENSEMBL/SYMBOL/ENTREZID
#if it is flybase or other ID you can convert the ID to ENSEMBL/SYMBOL/ENTREZID
#take only the ENSEMBL/SYMBOL/ENTREZID as input list so it will be used to search in database to find others 
#ENSEMBL/SYMBOL/ENTREZID but some of them are not found in Org.DB My MEL_0.05_2FC down from 6046 to 4853(20% unmapped)
#clusterprofiler need the ENSEMBL(MUST) and ENTREZID(VVVHelpful recommented)

df1 <- read.csv(file.choose())#male or female genes 
head(df1)
nrow(df1)


#for melanogaster only
names(df1)[1] <- "ENSEMBL"
names(df1)[2] <- "SYMBOL"
names(df1)[8] <- "log1FoldChange"


#No need of conversion if you have ENSEMBL as SAGD 
#you need entrezid ONLY for the kegg analysis 
#entrezid = ncbi id

#df1 <- subset(df1,select=c(ENSEMBL))
#df1$ENSEMBL
gene <- df1$ENSEMBL

keytypes(org.Dm.eg.db)

gene.df <- bitr(gene, fromType = "ENSEMBL", toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Dm.eg.db)
#gene.df <- bitr(gene, fromType = "ENSEMBL", toType = c("ENTREZID", "SYMBOL", "UNIPROT"), OrgDb = org.Dm.eg.db)

head(gene.df)
nrow(df1)
nrow(gene.df)

df1_Entrez <- gene.df
head(df1_Entrez)
head(df1)

#merge the main file with EntrezID for kegg analysis 
head(df1_Entrez)
head(df1)
nrow(df1_Entrez)
nrow(df1)

library(dplyr)
df3 <- left_join(df1,df1_Entrez) #matched based on commone name col There sould be one only #dplyr function
head(df3)
nrow(df3)

#need to remove NA
#this df3 Entrez has NA in the colume
library(tidyr)
df4 <- df3 %>% drop_na(ENSEMBL)
df4 <- df3 %>% drop_na(SYMBOL)
df4 <- df3 %>% drop_na(ENTREZID)
head(df4)
nrow(df4)
lost_gene <- nrow(df1)-nrow(df4)
lost_gene


###****
#DF4 HAS NO na IN eNTREZ ID AND READ FOR THE ENRICHMENT ANALYSIS 
###***

setwd("C:/Users/Mursalin/Desktop/SAGD Pub 1/SAGD Data sets_current/GO_Erichment/data for GSEA")
write.csv(gene.df, "Background_MEL_All_No_Filter_15946_Ens.Entz.Sym.csv", row.names = F)

#End OF THE PART 1 Male or Female DEGs clean up 

#WORK WITH BACKGROUND 

df2 <- read.csv(file.choose())#background
head(df2)
nrow(df2)


#for background
names(df2)[1] <- "ENSEMBL"
names(df2)[2] <- "SYMBOL"
names(df2)[8] <- "log1FoldChange"
head(df2)

#df2 is ready for the enrichGO













