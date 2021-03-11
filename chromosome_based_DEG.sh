#!/bin/bash
## To find the chromosome specific male and female biased gene expression number form the cleaned SAGD files
## Path is C:\Users\Mursalin\Desktop\SAGD Pub 1\SAGD Data sets_current\FC1-0.05
##

FILEID='Final_MEL_0.05_1.csv'

#awk -F',' 'FNR>1 && $7 ~ /^X/ && $8 ~/^-/{count++} END {print count+0}' $FILEID

#awk -F',' 'FNR>1 && $7 ~ /^X/ && $8!~/^-/{count++} END {print count+0}' $FILEID

array=( X Y 2L 2R 3L 3R 4)
for i in "${array[@]}"
do
echo $i
awk -F',' 'FNR>1 && $7 ~ /^'$i'/ && $8 ~/^-/{count++} END {print count+0}' $FILEID #negative log Female
awk -F',' 'FNR>1 && $7 ~ /^'$i'/ && $8!~/^-/{count++} END {print count+0}' $FILEID #positive log Male
done
