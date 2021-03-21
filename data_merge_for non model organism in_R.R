library(devtools)
library(ggplot2)
library(dplyr)


df1 <- read.csv(file.choose()) #the sim SAGD file
head(df1)
nrow(df1)
names(df1)[1] <- "ENSEMBL_DSim"
names(df1)[2] <- "SYMBOL_DSim"
names(df1)[8] <- "log1FoldChange"

df2 <- read.csv(file.choose()) #ensemble database file 
head(df2)
nrow(df2)
names(df2)[1] <- "ENSEMBL_DSim"
names(df2)[3] <- "ENSEMBL_DMel"
names(df2)[4] <- "SYMBOL_DMel"

#only selected col we need 
df2 <- subset(df2,select=c(ENSEMBL_DSim,ENSEMBL_DMel,SYMBOL_DMel))

#Only Unique gene in the file Ensemble metazoa file
df2_deduped <- df2[!duplicated(df2[c("ENSEMBL_DSim")]),] #removal of duplication based on col name 
#matched with ensemble metazoa count (13098 / 15546 Sim Genes have 13374/16349 Yak Genes
#homolog with Dmel)
head(df2_deduped)
nrow(df2_deduped)

df3 <- left_join(df1,df2_deduped) #matched based on commone name col There sould be one only #dplyr function
head(df3)
nrow(df3)
nrow(df1)
tail(df3)

#just to remove the names based on the col match it will not merge the data frames 
#df4 = df1[df2_deduped$ENSEMBL_DSim %in% df1$ENSEMBL_DSim,]
#head(df4)
#nrow(df3)


setwd("C:/Users/Mursalin/Desktop/SAGD Pub 1/SAGD Data sets_current/GO_Erichment/data for GSEA")
#write.csv(df3, "Final_Sim_P0.05_1FC_Male_Dmel_ID_added.csv", row.names = F)


#there many sim genes do not have Dmel orthologs It is better to remove those for futhrer use 
library(tidyr)
#tidyr has a new function drop_na
#based on col name 
df5 <- df3 %>% drop_na(ENSEMBL_DMel)
df5 <- df3 %>% drop_na(SYMBOL_DMel)
df5 <- df3 %>% drop_na(ENSEMBL_DMel,SYMBOL_DMel)
head(df5)
tail(df5)
nrow(df5)
lost_gene <- nrow(df1)-nrow(df5)
lost_gene
write.csv(df5, "Background_Final_Sim_with_Dmel_Ensemble_Symbol_removed_NA.csv", row.names = F)

#write.csv(df1_female, "Temp_MEL_P0.05_1FC_Female.csv", row.names = F)

#Only Unique gene in the file

# we want the log2 fold change 
original_gene_list <- df2$ENSEMBL_DSim

# name the vector
#names(original_gene_list) <- df$X
names(original_gene_list) <- df2$ENSEMBL_DMel

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
head(gene_list)
head(original_gene_list)
