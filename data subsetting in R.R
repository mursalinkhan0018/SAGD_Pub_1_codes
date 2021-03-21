#Data cleaning in R Subsetting 


library(devtools)
library(ggplot2)
df1 <- read.csv(file.choose())
head(df1)
nrow(df1)
names(df1)[1] <- "ENSEMBL"
names(df1)[2] <- "SYMBOL"
names(df1)[8] <- "log1FoldChange"


#only positive sign is male 
df1_male <- subset(df1, log1FoldChange >0 ,
              select=c(ENSEMBL:Padj))

#negative sign is making it female 
df1_female <- subset(df1, log1FoldChange <0 ,
                   select=c(ENSEMBL:Padj))

nrow(df1_male) + nrow(df1_female)

head(df1_male)
head(df1_female)


setwd("C:/Users/Mursalin/Desktop/SAGD Pub 1/SAGD Data sets_current/GO_Erichment/data for GSEA")
write.csv(df1_male, "Background_MEL_P0.05_1FC_Male.csv", row.names = F)
#write.csv(df1_female, "Temp_MEL_P0.05_1FC_Female.csv", row.names = F)

#make it poisitive 

df1_female_positive_log1FC <- df1_female[,8] * (-1)
write.csv(df1_female_positive_log1FC, "Temp_MEL_P0.05_1FC_Female_Pos_logFC.csv", row.names = F)
df2 <- read.csv(file.choose())
head(df2)
names(df2)[1] <- "log1FoldChangePos"

P0.05_1FC_Female_with_Pos_logFC <- cbind(df1_female, df2)

head(P0.05_1FC_Female_with_Pos_logFC)
nrow(P0.05_1FC_Female_with_Pos_logFC)
write.csv(P0.05_1FC_Female_with_Pos_logFC, "Final_MEL_P0.05_1FC_Fe_Pos_logFC.csv", row.names = F)

#take only a col data
df1_female <- subset(df1, log1FoldChange <0 ,
                     select=c(ENSEMBL:Padj))

nrow(df1_female)
df1_female_positive_log1FC_ENSEMBL <- df1_female[,1]
write.csv(P0.05_1FC_Female_with_Pos_logFC, "Final_MEL_P0.05_1FC_Fe_DAVID.csv", row.names = F)
