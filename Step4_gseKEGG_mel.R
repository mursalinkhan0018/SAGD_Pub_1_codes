#gseKEGG

library(ggplot2)
library(DOSE)
library(pathview)
library("AnnotationDbi")
library("org.Dm.eg.db")
library("ggnewscale")
library(GOSemSim)
library(enrichplot)
library(clusterProfiler)

##KEGG Gene Set Enrichment Analysis
#READ
#VVVVVI READ
#For KEGG pathway enrichment using the gseKEGG() function, we need to convert id types. 
#We can use the bitr function for this (included in clusterProfiler). 
#It is normal for this call to produce some messages / warnings.
#In the bitr function, the param fromType should be the same as keyType from the gseGO 
#function above (the annotation source).
#This param is used again in the next two steps: creating dedup_ids and df2.

#toType in the bitr function has to be one of the available options from keyTypes(org.Dm.eg.db) and must map to one of 'kegg',
#'ncbi-geneid', 'ncib-proteinid' or 'uniprot' because gseKEGG() only accepts one of these 4 options as it's keytype parameter.
#In the case of org.Dm.eg.db, none of those 4 types are available,
#but 'ENTREZID' are the same as ncbi-geneid for org.Dm.eg.db so we use this for toType.


#melanogater use only 

df1 <- read.csv(file.choose()) #male or female genes 
head(df1)
nrow(df1)

names(df1)[1] <- "ENSEMBL"
names(df1)[2] <- "SYMBOL"
names(df1)[8] <- "log1FoldChange"
names(df1)[10] <- "log1FoldChangePos" #female

# we want the log2 fold change in the input data 
#male
original_gene_list <- df1$log1FoldChange #log 1 FC or  2Fold changes
#female
original_gene_list <- df1$log1FoldChangePos #female

# make a vector from the data frame 
names(original_gene_list) <- df1$ENSEMBL

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (MUST for for GSEA analysis in clusterProfiler/KEGG)
gene_list = sort(gene_list, decreasing = TRUE)
head(gene_list)
head(original_gene_list)

#As our intial input, we use original_gene_list which we created above.

nrow(original_gene_list)

# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=org.Dm.eg.db)

#remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
head(dedup_ids)

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
###Make sure you are making newd2 for KEGGONLY
df2 = df1[df1$ENSEMBL %in% dedup_ids$ENSEMBL,] #VVI match and remove like awk
head(df2)
nrow(df2)
nrow(df1)
lostgene <- (nrow(df1)-nrow(df2))


# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID #add of the ncbi_id in the main data only matched lost some Like left.join
head(df2)
nrow(df2)


# Create a vector of the gene unuiverse NOT DATA Frame
#male
kegg_gene_list <- df2$log1FoldChange
#female 
kegg_gene_list <- df2$log1FoldChangePos


# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y


# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list) #VVI you need in future
head(kegg_gene_list)

# sort the list in decreasing order (required for GSEA clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

#Create gseKEGG object
#ncbi-geneid and ENTREZID are same :)

kegg_organism = "dme" #https://www.genome.jp/kegg/catalog/org_list.html #select dsi for simulans
kk <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")


kegg_summary <- data.frame(kk)
write.csv(kegg_summary, "gsekegg_Final_BP_MEL_FEMALE_0.05_1FC.csv", row.names = F)

#Now you can visulaized just like enrichGO and gseGO 

d <- godata('org.Dm.eg.db', ont="BP")
kk2 <- pairwise_termsim(kk, method = "JC", semData = d) #awesome 

dotplot(kk2, showCategory=30)

emapplot(kk2, showCategory=30)

gseaplot(kk2, by = "all", title = kk2$Description[1], geneSetID = 1)
#OR
gseaplot(kk2, geneSetID = 'dme01200')

#optional 

#Category Netplot

# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(kk2, categorySize="geneNum", foldChange=gene_list)

cnetplot(kk2, categorySize="pvalue", foldChange=gene_list)

#Ridgeplot
library("ggridges")
library(ggplot2)

ridgeplot(kk2) + labs(x = "enrichment distribution")

#save 600dpi plot without ggplot 
jpeg("Dotplot_Final_gseaKEGG_BP_MEL_FEMALE_0.05_1FC.jpeg", width = 16, height = 10, units = 'in', res = 600)
jpeg("Emapplot_Final_gseaKEGG_BP_MEL_FEMALE_0.05_1FC.jpeg", width = 20, height = 10, units = 'in', res = 600)
jpeg("Ridgeplot_Final_gseaKEGG_BP_MEL_FEMALE_0.05_1FC.jpeg", width = 16, height = 10, units = 'in', res = 600)



jpeg("GSEAplot_Final_gseaKEGG_BP_MEL_FEMALE_0.05_1FC_Des1.jpeg", width = 16, height = 10, units = 'in', res = 600)
jpeg("Cnetplot_Final_gseaKEGG_BP_MEL_FEMALE_0.05_1FC.jpeg", width = 20, height = 10, units = 'in', res = 600)
#any plot code will go here example #dotplot(ego, showCategory=50) 
dev.off()

