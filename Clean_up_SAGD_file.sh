#!/bin/bash
##No gene will be in the file less than 1 fold up or down regualted. New file is created for the clean data 
#To remove anything less than <=1 use ^-0. or ^0.(You just need this other fold change can download from SAGD :) )
#To remove anything less than <=2 use ^-1. or ^1.
#To remove anything less than <=3 use ^-2. or ^2.
#To remove anything less than <=4 use ^-3. or ^3.
#To remove anything less than <=8 use ^-7. or ^7.
#use one by one or create the loop 

## Define the tissue ID
FILEID="VIR" #no space in bash variable

awk -F "," '$8 ~ /^-0./ { next } { print }' ${FILEID}*.csv > temp1.csv #Negative log2(M/F) col 8 means Female biased 

awk -F "," '$8 ~ /^0./ { next } { print }' temp1.csv > Final_${FILEID}_0.05_1.csv #Positive log2(M/F) in log col 8 means Male biased 

