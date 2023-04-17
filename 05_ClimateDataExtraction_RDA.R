#Set directory
setwd("~/Desktop/Sarah/Salmon/COSEWIC/Manuscript_New_Figures/Climate_RDA/")

#Libraries
library(raster)
library(rgdal)
library(vegan)
library(ggrepel)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(rasterVis)
library(maps)
library(maptools)
library(mapdata)
library(rbioclim)
library(colorRamps)


#Get bioclim data - for present
bioclim <- recursive.getData(times="pres") # by default resolution 2.5m and the 19 bioclim variables. (but can be changed)
present_bioclim <- bioclim[["pres"]]

#Extract data for Bio 1-19
bio_present <- present_bioclim[[paste0("bio",c(1:19))]]

##Site locations - with lat/longs
all_loc <- read.csv("rivers_database.csv")

#Select specific columns needed
#Select DU name, site name, lat, and long
all_data1<- all_loc[,c(5, 7, 9,10)]

#Subset data for only Labrador rivers
all_data <- all_data1[which(all_data1$Proposed_DU=="DU02A" |
                             all_data1$Proposed_DU=="DU02B" |
                             all_data1$Proposed_DU=="DU02C" ), ]

#Create coordinate dataframe
coordinates <- cbind(all_data$Long, all_data$Lat)

#Get bioclim data for Labrador river/site coordinates
bioclim_all_sites_present <- raster::extract(bio_present, coordinates)

#Combine bioclim data with Site codes
bioclim_all_sites_present_info <- as.data.frame(cbind(as.character(all_data$Proposed_DU), bioclim_all_sites_present))

#Check for NAs -- may need to edit coordinates if any not on land
sum(is.na(bioclim_all_sites_present_info))

#Save data for all sites and Bioclim 1-19
write.table(bioclim_all_sites_present_info, "Enviro_data_Labrador.txt", quote=F, row.names = F, col.names = T, sep="\t")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

### Redundacy analysis (RDA) for Labrador data

#read in bioclimatic data for Labrador sites
results<- read.table("Enviro_data_Labrador.txt", header=T)

#Scale enviro data
scaleenv<-as.data.frame(apply(data.matrix(results[,2:20]), 2, function(x) scale(x)))


#Run RDA with bioclim/enviro data as response and genetic groups (or DUs) as contraining factor (V1)
salmon.rda <- rda(scaleenv~ results$V1, data=results)

#Results of RDA
summary(salmon.rda)
RsquareAdj(salmon.rda) #Note adjusted R2

#Quick plotting of RDA results
plot(salmon.rda, scaling=3)          # default is axes 1 and 2
plot(salmon.rda, choices = c(1, 3), scaling=3)  # axes 1 and 3



###Making a nice plot for RDA - include multiple steps here:

#I made these plots using ggplot by extracting information from the RDA (salmon.rda)

#Centroids of factors (DUs)
centroid_groups <- as.data.frame(scores(salmon.rda, display = 'cn', scaling=3))
centroid_groups$group<- rownames(centroid_groups)

#Fix names of groupings
centroid_groups$group[which(centroid_groups$group== "results$V1DU02A")] <- "NorthernLabrador"
centroid_groups$group[which(centroid_groups$group== "results$V1DU02B")] <- "LakeMelville"
centroid_groups$group[which(centroid_groups$group== "results$V1DU02C")] <- "SouthernLabrador"


#Extract point values for individuals
points_individuals <- as.data.frame(scores(salmon.rda, display="sites",  scaling=3))
points_individuals2 <- cbind(results[,1],points_individuals)
points_individuals2$DU <- points_individuals2$`results[, 1]`

#Extract loading values for bioclimatic data 
loading_clim <- as.data.frame(scores(salmon.rda, display="species",  scaling=3))
loading_clim$bioclim<-rownames(loading_clim)

#Plot all together
RDA_plot<- ggplot()+
  #This adds text for the bioclim variabels - uses geom_text_repel to avoid overlapping names
  geom_text_repel(data=loading_clim, aes(x=RDA1, y=RDA2, label=bioclim), col="skyblue4", size=4)+
#this adds arrows for the bioclim variabels (direction and length shows loadings on axes)
  geom_segment(data=loading_clim,  aes(x = 0, y = 0, xend = RDA1, yend = RDA2), col="skyblue2",
               arrow = arrow(length = unit(0.005, "npc")) )+
  #this adds the points for the individual salmon rivers - coloured by "DU"
  geom_point(data=points_individuals2, aes(x=RDA1,y= RDA2, bg=DU), pch=21, col="black", size=4)+
  
  #these are the colours for the groups (includes for pops and centroids) - (had to add colours twice??)
  scale_fill_manual(values=c("chartreuse3","#ffff33","lightpink2","#ffff33","chartreuse3","lightpink2"))+
  
  #This add the centroids for the factors (DUs)
  geom_point(data=centroid_groups, aes(x=RDA1,y= RDA2, fill=group), col="black", pch=24,  size=5)+
  
  #This adds the names for the factors (DUs)
  geom_text_repel(data=centroid_groups, aes(x=RDA1, y=RDA2, label=group), col="black", size=6)+

    #This adjusts the 'theme'/look of the plot
  theme_bw()+theme(panel.grid =element_blank(), axis.text = element_text(size=15), legend.position = "none")+

    #names for x and y axes (got proportions and R2 above)
  ylab("RDA2 (23.1%) ") + xlab("RDA1  (76.9%)")+
  
  geom_text(aes(x=1.2, y=1.2, label="R2adj=0.62\np=0.001"), cex=6)


#save plot
RDA_plot


#RDA Model significance
signif.full <- anova.cca(salmon.rda, parallel=getOption("mc.cores")) # default is permutation=999
signif.full
#By axis
signif.axis <- anova.cca(salmon.rda, by="axis", parallel=getOption("mc.cores"))
signif.axis

#Get loadings of bioclim variables
load.rda <- scores(salmon.rda, choices=c(1:2), display="species")  # Species scores for the first three constrained axes

#Save loadings of bioclim variables
loadings_variables<- as.data.frame(load.rda)
write.table(loadings_variables, "Loadings_RDA_Labrador.txt", quote = F, row.names = T, col.names = T, sep="\t")

#Check loadings to see which variables are loading highest on each axis
loadings_variables[order(abs(loadings_variables$RDA1), decreasing = T),]
loadings_variables[order(abs(loadings_variables$RDA2), decreasing = T),]


##############  Map of environmental data

#Note that BIO19 highest loading on RDA 1 and BIO10 highest loading on RDA 2

# plot your base graphics 
par(mfrow=c(1,2))

#extract bio19 data for plotting - downloaded raster of high res version from World clim
bio19 <- raster::raster("wc2.1_30s_bio_19.tif") 
#plot(bio19)

#get bio19 data for map (within lat/long specified)
r1_bio19_highrest <- crop(bio19, extent(-64, -54, 50.5,58))
#hist(r1_bio19_highrest)

#check values for bio19
values(r1_bio19_highrest)
#some high values are skewing plotting colours - so set these as NA as they are outside the geo region of interest anyway
values(r1_bio19_highrest)[values(r1_bio19_highrest) >300] = NA


#Create colour palette for plotting
colr <- colorRampPalette(brewer.pal(11, 'RdYlBu'))


#Plot map with bio19 and river locations - note often have to adjust window size in R Studio plotting
palette(c("#ffff33","chartreuse3","lightpink2")) #palette for salmon river points
# specify the color pallette
#Plot Bio19 raster data
plot(r1_bio19_highrest, xlab="Long", ylab="Lat",
     col= c(rep("midnightblue", 100),rep("#313695", 8000), rev(colr(10000)) ) ) 
#add map outline
map("worldHires", ylim=c(50.5, 58), xlim=c(-64,-54),fill=F, col="black",
    resolution=0, lwd=1, add=T);map.axes()
#add points for salmon rivers
points(x=all_data$Long,
       y=all_data$Lat, cex=1.5, 
       pch=21, col="black", bg=all_data$Proposed_DU)+title("bio19")


#extract bio10 data for plotting - downloaded raster of high res version from World clim
bio10<-raster::raster("wc2.1_30s_bio_10.tif") 

#plot(bio10)
r1_bio10_highrest <- crop(bio10, extent(-64, -54, 50.5,58))
#hist(r1_bio10_highrest)


#plot results for BIO10 on map
palette(c("#ffff33","chartreuse3","lightpink2")) # specify the color pallette for salmon river points
plot(r1_bio10_highrest, xlab="Long", ylab="Lat",
     col= c(rep("midnightblue", 10000), rep("#313695", 1000), rev(colr(10000)) ) ) 
#add map outline
map("worldHires", ylim=c(50.5, 58), xlim=c(-64,-54),fill=F, col="black",
    resolution=0, lwd=1, add=T);map.axes()
#add points for salmon rivers
points(x=all_data$Long,
       y=all_data$Lat, cex=1.5,
       pch=21, col="black", bg=all_data$Proposed_DU)+title("bio10")



