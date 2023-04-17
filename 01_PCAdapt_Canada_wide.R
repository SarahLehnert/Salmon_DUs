#Rscript to run PCAdapt on Canada-wide genomic data for Atlantic salmon

#set directory
setwd("~/Desktop/Sarah/Salmon/COSEWIC/Manuscript_New_Figures/")

#Read libraries
library(pcadapt)
library(ggplot2)
library(qqman)
library(qvalue)
library(ggrepel)
library(maps) # tool for maps
library(mapdata) # all your basemaps are here
library(ggman)
library(patchwork)

###Read in genotype data (plink format)
Canada_dat <- read.pcadapt("Combined_CAN_goodSNPs_allruns_badindremove_162K_MAF005_CanadaOnly_jan6_2023.bed", type="bed")

#Run pcadapt - use K=6 based on screeplot
Canada_dat_pca<- pcadapt(Canada_dat, K = 6, method =  "mahalanobis")

#see scree plot
plot(Canada_dat_pca, option = "screeplot")

#Get PCA scores for each sample
PCA_results <- Canada_dat_pca$scores
#Get variance explained by PC axes
pc_variance<-(Canada_dat_pca$singular.values)^2

#Get individual data (IDs and pop names) from .fam plink file
fam <- read.table("Combined_CAN_goodSNPs_allruns_badindremove_162K_MAF005_CanadaOnly_jan6_2023.fam", header=F)
#combine individual info with PCA results
PCA_results_CAN<- as.data.frame(cbind(fam, PCA_results))

#Update column names as needed
head(PCA_results_CAN)
colnames(PCA_results_CAN)=c("Pop", "ind", "x", "y", "z", "m", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6")

#Set pop as factor for plotting
PCA_results_CAN$Pop=as.factor(as.character(PCA_results_CAN$Pop))

#get addtional population info (lat/long and region of site locations)
pop_info<-read.csv("~/Desktop/Sarah/Salmon/COSEWIC/220K/List_Table_Samples_220K.csv", fill=T)
#Combine information with PCA results for plotting
combined_pop_PCA <- merge(PCA_results_CAN, pop_info[,c(1:3,6:8)], by=1)


#Plot PCA results with different regions highlighted by different colours
#Here highlight only Labrador as red - remaining locations are blue
ggplot()+geom_point(data=combined_pop_PCA, aes(PC1, PC2, fill=Region), col="black",pch=21, 
                               size=2.5)+
  scale_fill_manual(values=c("firebrick2", "midnightblue", "midnightblue", 
                             "midnightblue", "midnightblue", "midnightblue"))+theme_bw()


### Convert p-values to q-values for SNPs to identify outliers

#check number of SNPs
length(Canada_dat_pca$pvalues)

#Read in SNP name/locations
bim<-read.table("Combined_CAN_goodSNPs_allruns_badindremove_162K_MAF005_CanadaOnly_jan6_2023.bim", header=F)

#Combine pvalues with SNP info
bim_pvals <- as.data.frame(cbind(bim, Canada_dat_pca$pvalues ))

#Update column names as needed
head(bim_pvals)
colnames(bim_pvals)<- c("Chr", "SNP", "x", "Pos", "a", "b", "pval")

#convert p-values to q-values
convt_q<- qvalue(bim_pvals$pval)
convt_q$qvalues
bim_pvals$qvals<- convt_q$qvalues

#check minimum qvals
which.min(bim_pvals$qvals)

#check number of significant SNPs
nrow(bim_pvals[which(bim_pvals$qvals< 0.05),])

#Get list of outliers to save
outliers <- bim_pvals[which(bim_pvals$qvals< 0.05),]
nrow(outliers)
write.table(as.data.frame(outliers[,2]), file = "PCAdapt_significant_loci_CANADA_220k_K6.txt", col.names = F, row.names = F, sep="\t", quote = F)


#Run plink script to remove outliers to create neutral dataset
############## remove outleirs and run again:
system("cd ~/Desktop/Software/plink_mac_20200219/;
       ./plink --bfile ~/Desktop/Sarah/Salmon/COSEWIC/Manuscript_New_Figures/Combined_CAN_goodSNPs_allruns_badindremove_162K_MAF005_CanadaOnly_jan6_2023 --exclude ~/Desktop/Sarah/Salmon/COSEWIC/Manuscript_New_Figures/PCAdapt_significant_loci_CANADA_220k_K6.txt --make-bed --allow-extra-chr --chr-set 30 --out ~/Desktop/Sarah/Salmon/COSEWIC/Manuscript_New_Figures/Combined_CAN_goodSNPs_allruns_badindremove_162K_MAF005_CanadaOnly_jan6_2023_no_pca_outliers")

# repeat above script for PCAdapt to perform analysis on 'neutral datasets'
