#GSEGO analysis for melanogater 

#GSEGO analysis 
#REF #https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/

#melanogater use only 

df1 <- read.csv(file.choose()) #male or female genes 
head(df1)
nrow(df1)

names(df1)[1] <- "ENSEMBL"
names(df1)[2] <- "SYMBOL"
names(df1)[8] <- "log1FoldChange"
names(df1)[10] <- "log1FoldChangePos"

# we want the log2 fold change in the input data 
#male
original_gene_list <- df1$log1FoldChange #male #log 1 FC or  2Fold changes
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

#gene_list is ready to go into gseGO

#GSEA: Gene Set Enrichment
#Check which options are available with the keytypes command, for example 

keytypes(org.Dm.eg.db)

gse <- gseGO(geneList=gene_list, 
             ont ="BP", # or you can use ALL MF CC
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Dm.eg.db, 
             pAdjustMethod = "none")


cluster_summary <- data.frame(gse)
colnames(cluster_summary)
#setwd to save using session 
write.csv(cluster_summary, "Final_gseGO_BP_MEL_FEMALE_0.05_1FC.csv", row.names = F)


##Visualizing enrichGO clusterProfiler results 
#REF: #https://bioconductor.org/packages/devel/bioc/manuals/enrichplot/man/enrichplot.pdf
#Alhamdulliah

#dotplot_gene Ration in x

dotplot(gse, showCategory=30) #use count or gene ratio in x axis 

d <- godata('org.Dm.eg.db', ont="BP")
gse2 <- pairwise_termsim(gse, method = "JC", semData = d) #awesome 

emapplot_cluster(gse2,showCategory = 30)
emapplot(gse2, showCategory = 30) #Alhamdulliah

gseaplot(gse2, by = "all", title = gse$Description[1], geneSetID = 1)
gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)

# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(gse, categorySize="geneNum ", foldChange=gene_list, showCategory = 5,vertex.label.font=6)

#ridgeplot
install.packages("ggridges")
library("ggridges")
library(ggplot2)

ridgeplot(gse2) + labs(x = "enrichment distribution")


#save 600dpi plot without ggplot 
jpeg("Dotplot_Final_gseaGO_BP_MEL_FEMALE_0.05_1FC.jpeg", width = 16, height = 10, units = 'in', res = 600)
jpeg("Emapplot_Final_gseaGO_BP_MEL_FEMALE_0.05_1FC.jpeg", width = 20, height = 10, units = 'in', res = 600)
jpeg("Ridgeplot_Final_gseaGO_BP_MEL_FEMALE_0.05_1FC.jpeg", width = 16, height = 10, units = 'in', res = 600)



jpeg("GSEAplot_Final_gseaGO_BP_MEL_FEMALE_0.05_1FC_Des1.jpeg", width = 16, height = 10, units = 'in', res = 600)
jpeg("Cnetplot_Final_gseaGO_BP_MEL_FEMALE_0.05_1FC.jpeg", width = 20, height = 10, units = 'in', res = 600)
#any plot code will go here example #dotplot(ego, showCategory=50) 
dev.off()


