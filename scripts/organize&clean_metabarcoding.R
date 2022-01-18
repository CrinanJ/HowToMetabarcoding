##################################################
## Project: metabarcoding and networks
## Script purpose: function to organizing and cleaning metabarcoding data for statistical analysis 
## Date: 10/09/2021
## Author: Diogo F. Ferreira (ferreiradfa@gmail.com)
## Packages: dplyr and reshape2  
## Notes: to read everything more fluidly go to  Tools>Global Options>Code>check the box "Soft-wrap R source files"Apply)
###################################################

##OUTPUT: function returns a plot with cutline to remove samples based on control reads and data frame with 19 columns: 
  #prey: OTU number.
  #predator: sample number.
  #weight: number of reads.
  #otu_phylum: OTU phylum.
  #otu_class: OTU class.
  #otu_order: OTU order.
  #otu_family: OTU family.
  #otu_genus: OTU genus.
  #otu_species: OTU species.
  #year: year when the sample was collected.
  #season: season when the sample was collected.
  #landscape: landscape where the sample was collected.
  #farm: farm where the sample was collected.
  #species: predator species names.
  #species_code: predator species code (first letter of predator genus in uppercase with first three letters of specific epithet in lowercase).
  #animal: if predator is a bird or bat.
  #proportion: proportion of an OTU in a sample based on the number of reads.
  #freq_farm: presence/absence of an OTU in a farm (e.g. for zbj maximum is 16).

##INPUT: function to read, organize and clean metabarcoding data from pipeline by Rachel et al. Needs Metabarcoding data and 
#data: metabarcoding data with samples starting with a "X" and control with "Control". 
#samples: sample_list with species information for each sample in data (lab.nbr,year,season,landscape,farm,species,species_code,animal). See line 117.
#remove_samples: TRUE if you want to remove samples that have less reads then the controls. Using "Control" to identify control columns. See line 54.
#otus_clean: 0 to 5 if you want to exclude MOTUs with less than 0% to 5% of total number of reads per sample. Default is 1 for 1% (e.g. add 0.01 for 0.01%). Values above 5 return error. See line 144.
#keep_class: only keeping targeted taxonomic classes (default is c("Arachnida","Insecta")). If NULL keeps all classes. See line 157.
#remove_NAorders: TRUE if I want to remove MOTUs that are not identified until order. See line 168.
#remove_NAfamily: TRUE if I want to remove MOTUs that are not identified until family. See line 179.
#desired_species: species that I want to keep - uses species_code (first letter of genus in uppercase with first three letters of specific epithet in lowercase). See line 197.

final_metbar <- function(data = NA, sample_list = NA, remove_samples=F,otus_clean=1, keep_class=c("Arachnida","Insecta"),remove_NAorders=T,remove_NAfamily=F,desired_species=NULL){
 
  #######ORGANISING AND CLEANING METABARCODING DATA#####################
  metbarc <- zbj #metabarcoding data
  class(metbarc)
  #glimpse(metbarc)
  
  #Check the data - rows are otus and numbers are samples
  dim(metbarc)
  rownames(metbarc)
  colnames(metbarc)
  
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
  metbarc_clean1 <- subset(metbarc_clean, select = c(names(metbarc_clean[1:(which(grepl("Control" , names(metbarc)))[1]-1)])))#selects columns 1st and "ControlC1"
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
  
  
  ####Merging otus info 
  ###subsetting original metabarcoding data to get taxa info
  otus<-data.frame(otu=rownames(metbarc),subset(metbarc, select=c(phylum,class,order,family,genus,species)))
  head(otus)
  
  #merging order with links
  links_otu <- merge(links1, otus, by.x ="prey",by.y="otu",all.x = TRUE)
  head(links_otu)
  dim(links_otu)
  
  library(dplyr)
  links_order <- links_otu%>%
    select(c(predator,prey,weight,phylum,class,order,family,genus,species))%>%
    rename(otu_phylum=phylum,otu_class=class,otu_order=order,otu_family=family,otu_genus=genus,otu_species=species)
  
  head(links_order)
  dim(links_order)
 
  ####Merging metabarcoding data with species and habitat data
  ###importing database with feaces samples ID
  samples <<- sample_list 
  
  ###merging habitat data and species ID with metabarcoding data
  info <- subset(samples, select = c(lab.nbr,year,season,landscape,farm,species,species_code,animal))
  head(info)
  
  links_habitat <<- merge(links_order, info , by.x ="predator",by.y="lab.nbr",all.x = TRUE)
  head(links_habitat)
  dim(links_habitat)
  
  
  ###Excluding otus with less than 1% of total number of reads per sample
  ##Creating column with proportion of otus per individual
  links_propor <- data.frame()
  for (i in 1:length(unique(links_habitat$predator))){
    id <- links_habitat[links_habitat$predator==unique(links_habitat$predator)[i],]#select linls from a specific individual
    id$total_reads <- sum(id$weight)#calculates total number of reads for that individual
    id$proportion <- id$weight/id$total_reads#calculates proportions
    links_propor <- rbind(links_propor, id)
  } 
  head(links_propor)
  dim(links_propor)
  
    #if true only keeps rows with more than 1% of reads
  if(otus_clean==0){
    links_clean <- links_propor
  }else if(otus_clean>5){
    stop("Threshold to high (only values from 0 to 5 are accepted)", call. = FALSE) 
  } else {
    links_clean <- subset(links_propor, proportion >= (otus_clean/100))
  }
  dim(links_clean)
  head(links_clean)
  
  ##only want to keep relevant prey taxa - only want to keep Arthropoda:Arachnida and Insecta for now
  levels(as.factor(links_clean$otu_phylum))
  levels(as.factor(links_clean$otu_class))
  
  if(is.null(keep_class)){
    links_clean1 <- links_clean
  } else{
    links_clean1 <-links_clean%>%
      filter(otu_class %in% keep_class)#default c("Arachnida","Insecta") - this can be changed to include other classes
  }
  dim(links_clean1)
  head(links_clean1)
  
  ##checking order to decide if I should pull only the records that got a hit to ORDER
  if(remove_NAorders==T){
    sum(is.na(links_clean1$otu_order))
    levels(as.factor(links_clean1$otu_order))
    links_clean2 <- subset(links_clean1, otu_order !="NA")
  } else {
    links_clean2 <- links_clean1
  }
  dim(links_clean2)
  
  ##checking family to decide if I should pull only the records that got a hit to family
  if(remove_NAfamily==T){
    sum(is.na(links_clean2$otu_family))
    levels(as.factor(links_clean2$otu_family))
    links_clean3 <- subset(links_clean2, otu_family !="NA")
  } else {
    links_clean3 <- links_clean2
  }
  dim(links_clean3)
  
  ###only keeping the predator species that I want 
  if(is.null(desired_species)){
    links_filter <- links_clean3
  } else{
    links_filter <- links_clean3%>%
    filter(species%in%desired_species)#e.g. c("Hipposideros ruber")
  }
  head(links_filter)
  dim(links_filter)
  
  ##final dataset
  return(links_filter)
}

