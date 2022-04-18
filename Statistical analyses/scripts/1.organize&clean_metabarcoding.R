##################################################
## Project: How to metabarcode (Rachel et al 2022)
## Script purpose: function to organizing and cleaning metabarcoding data for statistical analysis 
## Date: 17/03/2022
## Author: Diogo F. Ferreira (ferreiradfa@gmail.com)
## Packages: dplyr and reshape2  
## Notes: to read everything more fluidly go to  Tools > Global Options> Code > check the box "Soft-wrap R source files" > Apply
###################################################

##OUTPUT: 
  #function returns a plot with cut-line to remove samples based on control reads 
  #data frame with 19 columns: 
    #predator: sample number
    #prey: ASV number
    #weight: number of reads per ASV and sample
    #asv_phylum: ASV phylum
    #asv_class: ASV class
    #asv_order: ASV order
    #asv_family: ASV family
    #asv_genus: ASV genus
    #asv_species: ASV species.
    #year: year when the sample was collected
    #season: season when the sample was collected
    #landscape: landscape where the sample was collected
    #farm: farm where the sample was collected
    #species: predator species name
    #species_code: predator species code (first letter of predator genus in uppercase with first three letters of specific epithet in lowercase).
    #animal: if predator is a bird or bat
    #total_reads: total number of reads per sample
    #proportion: proportion of an ASV in a sample based on the number of reads
    
##INPUT: function to read, organize and clean metabarcoding data from Rachel et al 2022. Needs Metabarcoding data (organised the same way) and information associated to samples 
  #data: metabarcoding data with samples coming 1st and then controls with "Control" in the name
  #samples: sample_list with species information for each sample in data (lab.nbr,year,season,landscape,farm,species,species_code,animal). See line 115
  #remove_samples: TRUE if you want to remove samples that have less reads then the controls. Using "Control" to identify control columns. See line 49
  #keep_class: only keeping targeted taxonomic classes (default is c("Arachnida","Insecta")). If NULL keeps all classes. See line 129
  #asvs_clean: 0 to 5 if you want to exclude ASVs with less than 0% to 5% of total number of reads per sample. Default is 1 for 1% (e.g. add 0.01 for 0.01%). Values above 5 return error. See line 141
  #remove_NAorders: TRUE if I want to remove ASVs that are not identified to order. See line 163
  #remove_NAfamily: TRUE if I want to remove ASVs that are not identified to family. See line 173
  #desired_species: predator species that I want to keep - uses species_code (first letter of genus in uppercase with first three letters of specific epithet in lowercase). See line 183

final_metbar <- function(data = NA, sample_list = NA, remove_samples=F,keep_class=c("Arachnida","Insecta"),asvs_clean=1,remove_NAorders=T,remove_NAfamily=F,desired_species=NULL){
 
  #######ORGANISING AND CLEANING METABARCODING DATA#####################
  metbarc <- data #metabarcoding data
 
  ####Cleaning data 
  ##excluding samples that have less reads then the controls - using "Control" to get controls and selecting samples as all columns between 1st column and 1st column with "Control"
    reads_samples <- colSums(subset(metbarc, select = c(names(metbarc[1:(which(grepl("Control" , names(metbarc)))[1]-1)]))))#selects columns between 1st and 1st column with "Control", and then gives total reads per  samples
    n_samples <- length(reads_samples)#number of samples
    reads_control <- colSums(subset(metbarc, select = c(grepl("Control" , names(metbarc)))))##grepl matches a regular expression to a target. Getting total reads per control
    reads_control#controls and its number of reads
    
    reads_remove <- reads_samples[reads_samples < max(reads_control)]#samples to be removed
    samples_removed <-names(reads_remove)
    n_removed <- length(samples_removed)#number of samples that were removed
    barplot(sort(reads_samples),ylim = c(0, 2*max(reads_control)),xlab=paste("total n of samples:",n_samples,";","samples to be removed:",n_removed),ylab = "number of reads")
    abline(h=max(reads_control),col="red")# you can choose a more proper value based on the plot if you think you're removing to much samples
    
    if(remove_samples==T){
      
    metbarc_clean <- subset(metbarc, select = setdiff(names(metbarc), samples_removed))
    
  } else {
    metbarc_clean <- metbarc
  }
  
  dim(metbarc_clean)
  
  ##ditching non-samples columns and the columns representing control samples
  metbarc_clean1 <- subset(metbarc_clean, select = c(names(metbarc_clean[1:(which(grepl("Control" , names(metbarc)))[1]-1)])))#selects columns until "ControlC1"
  dim(metbarc_clean1)
  
  
  ###Transform the bipartite matrices into a list (data-frame with links between predator and prey)
  library(reshape2)#for melt function
  metbarc.t <- data.frame(t(metbarc_clean1)) #transpose so I have individuals as id variable (rownames)
  rownames(metbarc.t)
  colnames(metbarc.t)
  metbarc.t2 <- cbind.data.frame(reference=row.names(metbarc.t),metbarc.t)#creating a reference to use in melt function
  links <- melt(metbarc.t2, na.rm = T)#changes organization of data - uses individuals as ID
  head(links)
  colnames(links) <- c("predator", "prey", "weight")
  links[,1]=as.character(paste(links[,1]))
  links[,2]=as.character(paste(links[,2]))
  head(links)
  dim(links)
  
  #Exclude all links with value 0 (zero) - need always to be removed due to melt function
  links1 <- subset(links, weight > 0)
  head(links1)
  dim(links1)
  
  
  ####Merging asvs info 
  ###subsetting original metabarcoding data to get taxa info
  asvs<-data.frame(asv=rownames(metbarc),subset(metbarc, select=c(phylum,class,order,family,genus,species)))
  head(asvs)
  
  #merging order with links
  links_asv <- merge(links1, asvs, by.x ="prey",by.y="asv",all.x = TRUE)
  head(links_asv)
  dim(links_asv)
  
  library(dplyr)
  links_order <- links_asv%>%
    select(c(predator,prey,weight,phylum,class,order,family,genus,species))%>%
    rename(asv_phylum=phylum,asv_class=class,asv_order=order,asv_family=family,asv_genus=genus,asv_species=species)
  
  head(links_order)
  dim(links_order)
 
  ####Merging metabarcoding data with species and habitat data
  ###importing database with feces sample ID
  samples <<- sample_list 
  
  ###merging habitat data and species ID with metabarcoding data
  info <- subset(samples, select = c(lab.nbr,year,season,landscape,farm,species,species_code,animal))
  head(info)
  
  links_habitat <<- merge(links_order, info , by.x ="predator",by.y="lab.nbr",all.x = TRUE)
  head(links_habitat)
  dim(links_habitat)
  
  ##only want to keep relevant prey taxa - only want to keep Arthropoda:Arachnida and Insecta for now
  levels(as.factor(links_habitat$asv_phylum))
  levels(as.factor(links_habitat$asv_class))
  
  if(is.null(keep_class)){
    links_clean <- links_habitat
  } else{
    links_clean <-links_habitat%>%
      filter(asv_class %in% keep_class)#default c("Arachnida","Insecta") - this can be changed to include other classes
  }
  dim(links_clean)
  head(links_clean)
  
  
  ###Excluding asvs with less than x% of total number of reads per sample
  ##Creating column with proportion of asvs per individual
  links_propor <- data.frame()
  for (i in 1:length(unique(links_clean$predator))){
    id <- links_clean[links_clean$predator==unique(links_clean$predator)[i],]#select links from a specific sample
    id$total_reads <- sum(id$weight)#calculates total number of reads for that sample
    id$proportion <- id$weight/id$total_reads#calculates proportions
    links_propor <- rbind(links_propor, id)
  } 
  head(links_propor)
  dim(links_propor)
  
    #if true only keeps rows with more than x% of reads
  if(asvs_clean==0){
    links_clean1 <- links_propor
  }else if(asvs_clean>5){
    stop("Threshold to high (only values from 0 to 5 are accepted)", call. = FALSE) 
  } else {
    links_clean1 <- subset(links_propor, proportion >= (asvs_clean/100))
  }
  dim(links_clean1)
  head(links_clean1)
  
  ##checking order to decide if I should pull only the records that got a hit to ORDER
  if(remove_NAorders==T){
    sum(is.na(links_clean1$asv_order))
    levels(as.factor(links_clean1$asv_order))
    links_clean2 <- subset(links_clean1, asv_order !="NA")
  } else {
    links_clean2 <- links_clean1
  }
  dim(links_clean2)
  
  ##checking family to decide if I should pull only the records that got a hit to family
  if(remove_NAfamily==T){
    sum(is.na(links_clean2$asv_family))
    levels(as.factor(links_clean2$asv_family))
    links_clean3 <- subset(links_clean2, asv_family !="NA")
  } else {
    links_clean3 <- links_clean2
  }
  dim(links_clean3)
  
  ###only keeping the predator species that I want 
  if(is.null(desired_species)){
    links_filter <- links_clean3
  } else{
    links_filter <- links_clean3%>%
    filter(species_code%in%desired_species)#e.g. c("Hipposideros ruber")
  }
  head(links_filter)
  dim(links_filter)
  
  ##final dataset
  return(links_filter)
}

