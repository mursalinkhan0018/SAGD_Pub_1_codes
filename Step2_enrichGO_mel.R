#enrichGO

#background and DEGs FOR enrichGO

allOE_genes <- as.character(df2$ENSEMBL) #make it string/character Background can have all genes in the sample
sigOE_genes <- as.character(df4$ENSEMBL) #make it string/character #only Padj 0.05 and Specidfic Fold Change>=1 6046 
head(allOE_genes)
head(sigOE_genes)

## Run GO enrichment analysis 
library(clusterProfiler)

ego <- enrichGO(gene = sigOE_genes, 
                universe = allOE_genes,
                keyType = "ENSEMBL",
                OrgDb = org.Dm.eg.db, 
                ont = "BP", # or you can use ALL MF CC
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

cluster_summary <- data.frame(ego)

write.csv(cluster_summary, "Final_enrichGO_BP_MEL_MALE_0.05_1FC.csv", row.names = F)


##Visualizing enrichGO clusterProfiler results 
#REF: #https://bioconductor.org/packages/devel/bioc/manuals/enrichplot/man/enrichplot.pdf
#Alhamdulliah

#dotplot_gene Ration in x
dotplot(ego, showCategory=50) #use count or gene ratio in x axis 

#using count data only in x axis 

barplot(ego, x= "GeneRatio",font.size = 12, title = "Barplot_Final_enrichGO_BP_MEL_MALE_0.05_1FC",
        label_format = 30, showCategory = 20) # x="Count"

#save 600dpi plot without ggplot 
jpeg("Dotplot_Final_enrichGO_BP_MEL_MALE_0.05_1FC.jpeg", width = 16, height = 10, units = 'in', res = 600)
jpeg("Barplot_Final_enrichGO_BP_MEL_MALE_0.05_1FC.jpeg", width = 16, height = 10, units = 'in', res = 600)
jpeg("Emapplot_Final_enrichGO_BP_MEL_MALE_0.05_1FC.jpeg", width = 30, height = 20, units = 'in', res = 600)
jpeg("Cnetplot_Final_enrichGO_BP_MEL_MALE_0.05_1FC.jpeg", width = 30, height = 20, units = 'in', res = 600)

#any plot code will go here 
#example #dotplot(ego, showCategory=50) 

dev.off()

#enrichment plot for nonmodel s

d <- godata('org.Dm.eg.db', ont="BP")
ego2 <- pairwise_termsim(ego, method = "JC", semData = d) #awesome 
emapplot(ego2, showCategory = 50)
emapplot_cluster(ego2, showCategory = 50)

### enrichment GO plot emapplot_alhamdulliah
#takes time use test file at first showcatefory 10

emapplot(ego, showCategory = 50)

##category netplot 
#Alhamdulliah
#use the test file for the work as it takes 4hours to run on background #showcategory1

OE_foldchanges <- df4$log1FoldChange
names(OE_foldchanges) <- df4$SYMBOL

#ego data frame is using from enrichment analysis 
#with ego2 it worked alhamdulliah
cnetplot(ego2, 
         categorySize="pvalue", 
         showCategory = 5, #upto 10 make be doable 5 is okay for SAGDPUB1
         foldChange=OE_foldchanges, 
         vertex.label.font=6)

#Each Species and Each Sex will have 4 plots 
#End of EnrichGO plots 
########################################################################################################################

#now we have enrichGo,gsego and gsekeggGO 
#play around for best figures 

heatplot(ego2, showCategory = 1, foldChange=OE_foldchanges) #forego created in cnetworkmap

