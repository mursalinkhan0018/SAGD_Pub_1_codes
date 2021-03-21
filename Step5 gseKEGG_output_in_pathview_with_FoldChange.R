#pathview in the kegg with gene expression data 

library(pathview)

#pathway ID from Kegg summary STEP4

#Final Pathview with all the gene expression 
load(url("http://genomedata.org/gen-viz-workshop/pathway_visualization/pathview_Data.RData"))
fc.kegg.sigmet.p.up[grepl("hsa03430", rownames(fc.kegg.sigmet.p.up), fixed=TRUE),]
head(fc.kegg.sigmet.p.up)
#this one is equivalent to our gseKEGG.summary output.csv file 
#we need just he pathway ID

#VVVVI

#you need .fc file formalt where it will be exact same as kegg_gene_list but with kegg_gene_list.fc
#kegg_gene_list.fc will have 
#example: ENTREZID(ncbiID actually and FoldChange) in a values (vector I think)

#head(tumor_v_normal_DE.fc)
#7105      64102       8813      57147      55732       2268 
#0.2319479 -0.3580517  0.7033971 -0.2258820  1.0475126  0.3837560 


#FROM STEP4

# sort the list in decreasing order (required for GSEA clusterProfiler)
kegg_gene_list.fc = sort(kegg_gene_list.fc, decreasing = TRUE) #.fc is must in the example format 
head(kegg_gene_list.fc)
head(kegg_gene_list)

pathview(gene.data=tumor_v_normal_DE.fc, species="hsa", pathway.id="hsa03430")
head(tumor_v_normal_DE.fc)

pathview(gene.data=kegg_gene_list.fc,species = "dme", pathway.id="dme03008") #pathwya.ID from gseaKEGG_summary.csv output 
