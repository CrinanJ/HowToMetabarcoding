##################################################
## Project: metabarcoding and networks
## Script purpose: function to organizing and cleaning metabarcoding data for statistical analysis 
## Date: 10/09/2021
## Author: Diogo F. Ferreira (ferreiradfa@gmail.com)
## Packages: dplyr, reshape2 and openxlsx  
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

##INPUT: function to read, organize and clean metabarcoding data from pipeline created by Andreanna. 
#data: metabarcoding data PLANTS 3,5 and 7 primer are missing the controls, so removel_samples should be set to FALSE to avoid error. FWH primer has a samples called X353Comb, so I changed it manually in the excel file to X353comb (if not X353Comb would be considered as a control).
#remove_samples: TRUE if I want to remove samples that have less reads then the controls. DOUBLE CHECK THIS WITH NEW DATA. Control names can be different between datasets. I'm using using "C" to identify control columns. See line 50.
#otus_clean: TRUE if I want to exclude OTUs with less than 1% of total number of reads per sample (criteria commonly use to clean metabarcoding data) See line 130.
#keep_class: only keeping targeted taxonomic classes (default is c("Arachnida","Insecta")). If NULL keeps all classes. See line 172.
#remove_NAorders: TRUE if I want to remove OTUs that are not identified until the order. See line 183.
#remove_NAfamily: TRUE if I want to remove OTUs that are not identified until the family. See line 194.
#desired_species: species that I want to keep - uses species_code (first letter of genus in uppercase with first three letters of specific epithet in lowercase). Birds data doesn't have a code for now. See line 212.

final_metbar <- function(data = NA,remove_samples=T,otus_clean=T,keep_class=c("Arachnida","Insecta"),remove_NAorders=T,remove_NAfamily=F,desired_species=NULL){
 
  library(dplyr)#for glimpse and to organize c
  ###############################ORGANISING AND CLEANING METABARCODING DATA##############################
  metbarc <- zbj #choose desired primer
  class(metbarc)
  #glimpse(metbarc)
  
  #Check the data - rows are otus and X[number] are individuals
  dim(metbarc)
  rownames(metbarc)
  colnames(metbarc)
  
  ####Cleaning data 
  ##excluding samples that have less reads then the controls - I think the control names is different between datasets but using "C" should be okay - IN FWH DOESN'T WORK DUE TO X353Comb - DOUBLE CHECK THIS WITH NEW DATA
  if(remove_samples==T){
    reads_samples <- colSums(subset(metbarc, select = c(grepl("X" , names(metbarc)))))#reads total per  samples
    n_samples <- length(reads_samples)#number of samples
    reads_control <- colSums(subset(metbarc, select = c(grepl("C" , names(metbarc)))))#read total per control
    reads_control#controls and its number of reads
    
    reads_remove <- reads_samples[reads_samples < max(reads_control)]#samples to be removed
    samples_removed <-names(reads_remove)
    n_removed <- length(samples_removed)#number of samples that were removed
    barplot(sort(reads_samples),ylim = c(0, max(reads_control)),xlab=paste("total samples:",n_samples,";","samples removed:",n_removed),ylab = "number of reads")#I can choose a more proper value based on the plot if I think I'm removing to much samples - zbj has a control with 743 reads (it's a lot)
    
    metbarc_clean <- subset(metbarc, select = setdiff(names(metbarc), samples_removed)) #had to remove samples for all datasets (see "checking species to keep.r")
  } else {
    metbarc_clean <- metbarc
  }
  
  dim(metbarc_clean)
  
  ##ditching non-taxa columns and the columns representing control samples
  metbarc_clean1 <- subset(metbarc_clean, select = c(grepl("X" , names(metbarc_clean))))#grepl matches a regular expression to a target (X for only individuals)
  dim(metbarc_clean1)
  
  
  ###Transform the bipartite matrices into a edge list (data.frame with links between predator and prey)
  library(reshape2)#for melt function
  metbarc.t <- data.frame(t(metbarc_clean1)) #transpose so I have individuals as id variable (rownames)
  rownames(metbarc.t)
  colnames(metbarc.t)
  metbarc.t2 <- cbind.data.frame(reference=row.names(metbarc.t),metbarc.t)#creating a reference to use in melt function
  links <- melt(metbarc.t2, na.rm = T)#changes organization of data - uses individuals as ID
  colnames(links) <- c("predator", "prey", "weight")
  links[,1]=as.character(paste(links[,1]))#
  links[,2]=as.character(paste(links[,2]))
  head(links)
  dim(links)
  
  #Exclude all edges with value 0 (zero) - need always to be removed due to melt function
  links1 <- subset(links, weight > 0)
  head(links1)
  dim(links1)
  
  
  ####Merging otus info 
  ###subsetting original metabarcoding data to get taxa info
  otus<-data.frame(otu=rownames(metbarc),subset(metbarc, select=c(phylum,class,order,family,genus,species,sacc)))
  head(otus)
  
  #merging order with links
  links_otu <- merge(links1, otus, by.x ="prey",by.y="otu",all.x = TRUE)
  head(links_otu)
  dim(links_otu)
  
  links_order <- links_otu%>%
    select (predator,prey,weight,phylum,class,order,family,genus,species)%>%
    rename(otu_phylum=phylum,otu_class=class,otu_order=order,otu_family=family,otu_genus=genus,otu_species=species)
  
  head(links_order)
  dim(links_order)
 
  ####Merging metabarcoding data with species and habitat data
  ###importing database with feaces samples ID
  library(openxlsx)#read excel and sheet
  samples <<- read.xlsx("data/metabarcoding/feces_sample_database.xlsx", sheet = "samples")
  
  ##when imported to R metabarcoding id samples in metabarcoding data have an X before number
  samples$lab.nbr <- paste0("X",samples$lab.nbr)#add X to lab numbers
  head(samples)
  
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
    id <- links_habitat[links_habitat$predator==unique(links_habitat$predator)[i],]#select edges from a specific individual
    id$total <- sum(id$weight)#calculates total number of reads for that individual
    id$proportion <- id$weight/id$total#calculates proportions
    links_propor <- rbind(links_propor, id)
  } 
  head(links_propor)
  dim(links_propor)
  
    #if true only keeps rows with more than 1% of reads and otus that appear in at least half of the farms
  if(otus_clean==T){
    links_clean <- subset(links_propor, proportion >= 0.01)
  } else {
    links_clean <- links_propor
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
      filter(otu_class %in% keep_class) #c("Arachnida","Insecta") - this can be changed to include other classes
  }
  dim(links_clean1)
  head(links_clean1)
  
  ##checking order to decide if I should pull only the records that got a hit to ORDER
  #Too many NA (need to fix this in the future), but for now will use it to have a smaller dataset
  if(remove_NAorders==T){
    sum(is.na(links_clean1$otu_order))
    levels(as.factor(links_clean1$otu_order))
    links_clean2 <- subset(links_clean1, otu_order !="NA")
  } else {
    links_clean2 <- links_clean1
  }
  dim(links_clean2)
  
  ##checking family to decide if I should pull only the records that got a hit to family
  #Too many NA (need to fix this in the future), but for now will use it to have a smaller dataset
  if(remove_NAfamily==T){
    sum(is.na(links_clean2$otu_family))
    levels(as.factor(links_clean2$otu_family))
    links_clean3 <- subset(links_clean2, otu_family !="NA")
  } else {
    links_clean3 <- links_clean2
  }
  dim(links_clean3)
  
  ###see if there is missing data - X146 doesn't have info but it was removed due to low nbr of reads
  #links_clean3[is.na(links_clean3$species) > 0,]
  
  #remove individuals that don't have data for now - need to add the missing data
  #links_clean3 <- links_clean3[links_clean3$predator!="X146",]
  #dim(links_clean3)
  #head(links_clean3)
  
  ###only keeping the predator species that I want 
  if(is.null(desired_species)){
    links_filter <- links_clean3
  } else{
    links_filter <- links_clean3%>%
    filter(species_code%in%desired_species)
  }
  head(links_filter)
  dim(links_filter)
  
  ##final dataset
  return(links_filter)
}

