#script for running parallel structure with 101 microsatellites for Labrador
#same script used to run 96-SNP dataset in parallel structure

#Run script in terminal
library(parallel)
library(ParallelStructure)
library(genepopedit)

#Set directory
setwd("~/Desktop/Sarah/Salmon/COSEWIC/Micro/DUs/DU2_Labrador/Panel_101micr/")

#Get pop/individual info to create a job file - use genepop file to get info
all_inf<-genepop_detective("Lab_AtlSalmon_Micros_noPR_editgenos.gen", "All")

#Create structure file from genepop using genepopedit
genepop_structure(genepop = "Lab_AtlSalmon_Micros_noPR_editgenos.gen",locusnames = F,
                  path = "Lab_AtlSalmon_Micros_noPR_editgenos_structure.txt")

#set new directory to output structure
setwd("~/Desktop/Sarah/Salmon/COSEWIC/Micro/DUs/DU2_Labrador/Panel_101micr/Structure/")

#create joblist file (see function script at end of this script)
#run based on number of pops and K
JobList(nPops = 34, k=10, dir="out")

#Path to structure
STR_path="~/Desktop/Software/structure/"


#Run parallel structure
#info for individuals, loci, etc is from genepop_detective output

parallel_structure(structure_path=STR_path, 
                   joblist="out34_10_500000_100000_3.txt",
                   n_cpu=30, 
                   infile="Lab_AtlSalmon_Micros_noPR_editgenos_structure.txt",  onerowperind = 0,
                   outpath="Results/", numinds=length(all_inf$Inds), numloci=length(all_inf$Loci), printqhat=1,
                   plot_output=1, popdata=1, label=1, markernames=0, usepopinfo=0)

##after structure is complete

##### Zip files for structure harvester
system("cd ~/Desktop/Sarah/Salmon/COSEWIC/96SNP/DUs/DU2_Labrador/STRUCTURE/
       zip -vr all.zip results -x *.pdf")

##### Zip files for clumpak
setwd("~/Desktop/Sarah/Salmon/COSEWIC/96SNP/DUs/DU7_QuebecEastNShore/STRUCTURE/Results/")
currentdir <- getwd()
dir.create("ffiles")

newdir <- "~/Desktop/Sarah/Salmon/COSEWIC/96SNP/DUs/DU7_QuebecEastNShore/STRUCTURE/Results/ffiles"

files <- list.files(path = currentdir,  pattern = "_f", full.names = TRUE)

files_new <- gsub(dirname(files[1]), newdir, files)

for (i in 1:length(files)) {
  
  file.copy(files[i], files_new[i])
}

system("cd ~/Desktop/Sarah/Salmon/COSEWIC/96SNP/DUs/DU7_QuebecEastNShore/STRUCTURE/Results/
  zip -vr ffiles.zip ffiles -x *.pdf")

# delete the directory -- must add recursive = TRUE
unlink("~/Desktop/Sarah/Salmon/COSEWIC/96SNP/DUs/DU7_QuebecEastNShore/STRUCTURE/Results/ffiles", recursive = TRUE)
#Run in clumpak



######################

#JobList function from Ryan Stanley for parallel structure

JobList <- function(nPops,k,nsim=500000,burnin=100000,reps=3,dir=NULL){
  
  #This function will create a jobslist which can be used by paralell structure
  
  #   npops is the number of populations in your data
  #   k is a number or vector of potentail clustering that STRUCTURE will look for
  #     if k = 5 then it will look for 1 through 5 clusters
  #     if k = c(1,3,4) it will look for 1, 3 and 4 clusters
  #   nsim is the number of simulations (defaults to 500000)
  #   burnin is the burnin desired for MCMC (defaults to 100000)
  #   reps is the number of replicate clusterings for each level of k (defaults to 3)
  #   dir is the directory you where you want the job list to be created. By default it will just return
  #     to R environment
  
  
  options(scipen = 999) # this will stop R from displaying the burning and nsim values in scientic notation
  nsim=as.numeric(nsim)
  
  #Make a population vector
  if(length(k)==1)
  {
    kvec <- rep(1:k,each=reps)
  }
  
  if(length(k)>1)
  {
    kvec <- rep(k,each=reps)
  }
  
  #Population data
  temp1 <- rep(paste(as.character(1:nPops), collapse=","))
  PopVec <- as.vector(do.call("rbind", replicate(length(kvec),temp1, simplify = FALSE)))
  jobs <- paste("T",1:length(kvec),sep="")
  burn <- rep(burnin,length(kvec))
  sim <- rep(nsim,length(kvec))
  
  temp2 <- as.data.frame(cbind(jobs,PopVec,kvec,burn,sim))
  
  output <- do.call(paste,c(temp2, sep=" "))
  
  if(length(dir)>0)
  {
    write.table(output,paste(dir,nPops,"_",max(k),"_",nsim,"_",burnin,"_",reps,".txt",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE)
  }
  
  if(length(dir)<=0)
  {
    return(output)
  }
  
  
}
