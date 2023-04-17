#Script for running PCAdapt and Manhattan plot

#Set directory
setwd("~/Desktop/Sarah/Salmon/COSEWIC/220K/DUs/DU2_Labrador/Files/noPR/")

#Libraries
library(pcadapt)
library(ggplot2)
library(qqman)
library(qvalue)
library(ggrepel)

#Read in site names and location for all 
site_info<- read.csv("~/Desktop/Sarah/Salmon/CIGENE_data/Updated_Dataset/Final/Canada/All_siteLocations_NorthAmerican_UPDATED.csv", header=T)

#Read in plink genotype file for pcadapt
pca_adapt<- read.pcadapt("DU2_Lab_220K_nobadloc_maf005_noPR.bed", type="bed")
#run pcadapt
pca_adapt_out<- pcadapt(pca_adapt, K = 2, method =  "mahalanobis")
#check scree plot
plot(pca_adapt_out, option = "screeplot")

#Get variance explained by PC axes
pc_variance<-(pca_adapt_out$singular.values)^2

#Get individual scores
PCA_results <- pca_adapt_out$scores

#Read in plink .fam file to combine with pcadapt results
fam<-read.table("DU2_Lab_220K_nobadloc_maf005_noPR.fam", header=F)

#Combine PCAdapt results and individual info
PCA_results_DU12<- as.data.frame(cbind(fam, PCA_results))

#Check headings and update column names
head(PCA_results_DU12)
colnames(PCA_results_DU12)=c("Pop", "ind", "x", "y", "z", "m", "PC1", "PC2")

#Set pops as factor and add column for 'region' to specify the three DUs in labrador
PCA_results_DU12$Pop=as.factor(as.character(PCA_results_DU12$Pop))
PCA_results_DU12$Region=rep("new")

#Update which sites belong to which regions for plotting
PCA_results_DU12$Region[which(PCA_results_DU12$Pop=="MU" |
                                PCA_results_DU12$Pop=="SK" |
                                PCA_results_DU12$Pop=="CB" |
                                PCA_results_DU12$Pop=="KE" |
                                PCA_results_DU12$Pop=="TR" |
                                PCA_results_DU12$Pop=="CK" |
                                PCA_results_DU12$Pop=="RW" |
                                PCA_results_DU12$Pop=="SR" |
                                PCA_results_DU12$Pop=="MB" |
                                PCA_results_DU12$Pop=="CL" |
                                PCA_results_DU12$Pop=="CR" )]<- "LakeMelville"

PCA_results_DU12$Region[which(PCA_results_DU12$Pop=="HU" |
                                PCA_results_DU12$Pop=="ENG" |
                                PCA_results_DU12$Pop=="BIG" )]<- "NorthernLabrador"

PCA_results_DU12$Region[which(PCA_results_DU12$Pop=="PA" |
                                PCA_results_DU12$Pop=="SH" |
                                PCA_results_DU12$Pop=="LL" |
                                PCA_results_DU12$Pop=="CHR" |
                                PCA_results_DU12$Pop=="FO" |
                                PCA_results_DU12$Pop=="EA" )]<- "SouthernLabrador"


#Get mean PC values for each population on PC1 and PC2 - to use for plotting
PC1_mean<-aggregate(PCA_results_DU12$PC1~PCA_results_DU12$Pop, FUN = "mean")
PC2_mean<-aggregate(PCA_results_DU12$PC2~PCA_results_DU12$Pop, FUN = "mean")
aggregate_means<-as.data.frame(cbind(PC1_mean, PC2_mean[,2]))
colnames(aggregate_means)=c("Pop", "PC1_mean", "PC2_mean")

#PCA plot of Labrador genomic dataset
ggplot()+geom_point(data=PCA_results_DU12, aes(PC1, PC2, bg=Region), 
                   col="black", pch=21, size=3)+theme_bw()+
  scale_fill_manual(values=c("#ffff33","chartreuse3","lightpink2"))+
  xlab(paste0("PC 1 (",round(pc_variance[1],3)*100, "%)" )) +
  ylab(paste0("PC 2 (",round(pc_variance[2],3)*100, "%)" )) +
 geom_text_repel(data=subset(aggregate_means, PC1_mean> 0), aes(PC1_mean, PC2_mean, label=Pop),
                 nudge_x       = 0.01 + subset(aggregate_means, PC1_mean> 0)$PC1_mean,
                 segment.size  = 0.2,
                 segment.color = "grey50",
                 direction     = "y",
                 hjust         = 1,
                 col="black") +
  geom_text_repel(data=subset(aggregate_means, PC1_mean< -0.05), aes(PC1_mean, PC2_mean, label=Pop),
                  nudge_x      = 0.002 + subset(aggregate_means, PC1_mean< -0.05)$PC1_mean,
                  nudge_y     = 0.002+ subset(aggregate_means, PC1_mean< -0.05)$PC2_mean,
                    segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "y",
                  hjust         = 0.5,
                  col="black") +
  geom_text_repel(data=subset(aggregate_means, PC1_mean > -0.05 & PC1_mean < 0), aes(PC1_mean, PC2_mean, label=Pop),
                  nudge_x      = 0.01 - subset(aggregate_means, PC1_mean > -0.05 & PC1_mean < 0)$PC1_mean,
                  nudge_y     = 0.025+ subset(aggregate_means, PC1_mean > -0.05 & PC1_mean < 0)$PC2_mean,
                  segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "y",
                  hjust         = 1,
                  col="black") +
  theme(axis.title = element_text(size=15), axis.text=element_text(size=13), legend.text = element_text(size=13))
  NULL

###################################################
  
###Manhattan plotting of results

#Check number of SNPs
length(pca_adapt_out$pvalues)

#Read in SNP data file (.bim from plink)
bim<-read.table("DU2_Lab_220K_nobadloc_maf005_noPR.bim", header=F)
#combine SNP info with p-values from pcadapt
bim_pvals <- as.data.frame(cbind(bim, pca_adapt_out$pvalues ))

#update column names as needed
head(bim_pvals)
colnames(bim_pvals)<- c("Chr", "SNP", "x", "Pos", "a", "b", "pval")

#Convert p-values to qvalues
convt_q<- qvalue(bim_pvals$pval)
convt_q$qvalues
bim_pvals$qvals<- convt_q$qvalues

#Check minimu q-val for plotting
which.min(bim_pvals$qvals)
bim_pvals[order(bim_pvals$qvals, decreasing = F),]

#Remove any values that are NA as won't plot - shouldn't apply here
bim_pvals<- bim_pvals[!is.na(bim_pvals$qvals),]

#Remove any Chr that are NAs - shouldn't apply to any here
plot_dat<- as.data.frame(bim_pvals[!is.na(bim_pvals$Chr),])
colnames(plot_dat)

#Plot Manhattan using qqman
qqman::manhattan(x = plot_dat, cex=1,  col=c("gray50", "dodgerblue4"),
          chr = "Chr", bp = "Pos", snp = "SNP", ylab="-log10(qvalue)",
          p = "qvals", suggestiveline = F, genomewideline = -log10(0.05) )


#Check number of significant loci (outliers q<0.05)
nrow(bim_pvals[which(bim_pvals$qvals< 0.05),])
#314 significant loci

#Check which chromosomes they fall to
table(bim_pvals$Chr[which(bim_pvals$qvals< 0.05)])


##Save files with all loci and their qvalues, as well as a data file of all outlier loci
#these will be used for Gene Ontology analyses

outliers<-bim_pvals[which(bim_pvals$qvals< 0.05),]

write.table(outliers, file = "PCAdapt_significant_loci_du2_220k.txt", col.names = F, row.names = F, sep="\t", quote = F)
write.table(bim_pvals, file = "PCAdapt_allloci_du2_220k.txt", col.names = F, row.names = F, sep="\t", quote = F)


