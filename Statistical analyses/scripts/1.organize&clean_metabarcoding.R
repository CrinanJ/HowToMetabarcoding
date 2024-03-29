##################################################
## Project: Vogelwarte metabarcoding
## Script purpose: function to organizing and cleaning metabarcoding data for statistical analysis 
## Date: 17 August 2023
## Author: Diogo F. Ferreira (ferreiradfa@gmail.com), Crinan Jarrett (crinan.jarrett@vogelwarte.ch) & Lara Gross (lara.gross@vogelwarte.ch)
## Packages: dplyr and reshape2  
## Notes: to read everything more fluidly go to: Tools>Global Options>Code>check the box "Soft-wrap R source files">Apply)
###################################################

## OUTPUT: function a data frame with the following columns (if remove_samples=T, it also returns a plot with cutline to remove samples based on control reads): 
  #prey:            ASV number.
  #predator:        sample number.
  #weight:          number of reads.
  #asv_phylum:      ASV phylum.
  #asv_class:       ASV class.
  #asv_order:       ASV order.
  #asv_family:      ASV family.
  #asv_genus:       ASV genus.
  #asv_species:     ASV species.
  #site:            site where the sample was collected.
  #species:         predator species names.
  #total_reads:     total number of reads for that individual (sample) - AFTER cleaning.
  #total_asvs:      total number of ASVs for that individual (sample) - AFTER cleaning.
  #proportion:      proportion of an ASV in a sample based on the number of reads - AFTER cleaning.
  #reads_hit:       number reads for this taxon hit in a sample.
  #reads_hit_proportion: proportion of number reads for this taxon hit to all reads in a sample.
  
## INPUT: function to read, organize and clean metabarcoding data from pipeline by Andreanna Welch et al. 
  #data:            metabarcoding data with following format: samples starting with a letter prefix and control starting with "Control" (or other name). ASVs need to be read in as row names (e.g., read.csv("...", row.names="ASV")).
  #samples:         sample_list with species information for each sample in data (lab.nbr, animal and species COMPULSORY, can add other variables like site, age and sex OPTIONALLY). See line 188.
  #prefix_control:  in the data sheet, the prefix used to indicate control samples
  #remove_samples:  TRUE if you want to remove samples that have less reads then the controls. Using "Control" to identify control columns. See line 101.
  #remove_control_ASVs: TRUE if you want to remove ASVs that appear in control above cut-off (default: 1% of full run), FALSE if you want to keep all ASVs. See line 123.
  #keep_class:      only keeping targeted taxonomic classes (default is c("Arachnida","Insecta")). If NULL keeps all classes. See line 196.
  #hits_clean:      0 to 5; if you want to exclude taxon hits that make up less than 0% to 5% of total number of reads per sample. Default is 1 for 1% (e.g., add 0.01 for 0.01%). Values above 5 return error. See line 210.
  #asvs_clean:      0 to 5; if you want to exclude ASVs with less than 0% to 5% of total number of reads per sample. Default is 1 for 1% (e.g., add 0.01 for 0.01%). Values above 5 return error. See line 250.
  #control_ASVs_clean: percentage; if you want to remove ASVs that appear in control above the chosen cut-off (default: 1% of total reads of this ASV). Recommended. See line 276.
  #remove_NAorders: TRUE if I want to remove ASVs that are not identified until order. See line 306.
  #remove_NAfamily: TRUE if I want to remove ASVs that are not identified until family. See line 318.
  #desired_species: species that I want to keep - uses species name. See line 330.

## ! ##
## ! ## If the input data contains data about several metabarcoding runs (each with one or more controls), this function should be applied separately to each run.
## ! ##


final_metbar <- function(data = NA, sample_list = NA, prefix_control = "Control", 
                         remove_samples=F, remove_control_ASVs=F, keep_class=c("Arachnida", "Insecta"), hits_clean=0, asvs_clean=1,
                         control_ASVs_clean=1, remove_NAorders=F, remove_NAfamily=F, desired_species=NULL){
  # Packages needed --------------------------
  library(dplyr)
  library(reshape2) #for melt function
  library(tidyr) #for replace_na function
  options(scipen = 999) # turns off scientific notation, good for step "control_ASVs_clean"
  options(error=NULL) #I needed this, otherwise it opens an error browser when an error is thrown
  options(dplyr.summarise.inform = FALSE) #Turn of the warnings when grouping variables for dplyr::summarise()
  
  # Organising & cleaning metabarcoding data --------------------------
  metbarc <- data #metabarcoding data
  
  ## If these columns are not present in sample_list, add them (for the sake of the code). They will be removed at the end if empty.
  if(!("site" %in% names(sample_list))) {
    sample_list$site <- NA
  }
  if(!("season" %in% names(sample_list))) {
    sample_list$season <- NA
  }
  if(!("year" %in% names(sample_list))) {
    sample_list$year <- NA
  }
  if(!("landscape" %in% names(sample_list))) {
    sample_list$landscape <- NA
  }
  if(!("collection_date" %in% names(sample_list))) {
    sample_list$collection_date <- NA
  }
  if(!("age" %in% names(sample_list))) {
    sample_list$age <- NA
  }
  if(!("sex" %in% names(sample_list))) {
    sample_list$sex <- NA
  }
  if(!("long" %in% names(sample_list))) {
    sample_list$long <- NA
  }
  if(!("lat" %in% names(sample_list))) {
    sample_list$lat <- NA
  }
  
  ## Check the data - rows are ASVs, numbers are samples
  class(metbarc)
  dim(metbarc)
  rownames(metbarc)
  colnames(metbarc)
  
  ## Identify sample prefix used in dataset
  prefix <- as.character(gsub("[^a-zA-Z]", "", names(metbarc)[2]))

  ## Excluding samples that have less reads than the controls - using "Control" to get controls and selecting samples as all columns between 1st column and 1st column with "Control"
  reads_samples <- colSums(dplyr::select(metbarc, starts_with(prefix))) #selects sample columns, and then gives total reads per samples
  n_samples <- length(reads_samples) #number of samples
  if('try-error' %in% class(try(colSums(dplyr::select(metbarc, starts_with(prefix_control))))) |
     length(colSums(dplyr::select(metbarc, starts_with(prefix_control)))) == 0){
    stop("Prefix of control column needs to be specified (e.g., prefix_control='EC' for control column EC1)", call. = FALSE)
  }
  reads_control <- colSums(dplyr::select(metbarc, starts_with(prefix_control))) #getting total reads per control
  reads_control #number of reads for each control
  reads_total <- sum(reads_samples) + reads_control
  
  
  # REMOVE_SAMPLES --------------------------
  ### If you want to remove samples that have less reads then the controls. Using "Control" to identify control columns.
  if(remove_samples == T){
    ## Which samples should be removed?
    reads_remove <- reads_samples[reads_samples < max(reads_control)] #samples to be removed
    samples_removed <-names(reads_remove)
    n_removed <- length(samples_removed) #number of samples that were removed
    
    ## Plot of samples that will be removed
    barplot(sort(reads_samples), ylim = c(0, 2*max(reads_control)), 
            xlab=paste("total n of samples:", n_samples,";","samples to be removed:", n_removed), 
            ylab = "number of reads")
    abline(h = max(reads_control), col="red") #you can choose a more proper value based on the plot if you think you're removing to much samples
    
    ## Remove
    metbarc_clean <- subset(metbarc, select = setdiff(names(metbarc), samples_removed))
  } else {
    metbarc_clean <- metbarc
  }
  dim(metbarc_clean)
  
  
  # REMOVE_CONTROL_ASVS --------------------------
  ### Remove ASVs from dataset when they appear in controls at >1% of total run reads. Total run reads are reads from: all ASVs, all samples, all controls.
  if(remove_control_ASVs == T){
    control_props <- dplyr::select(metbarc_clean, starts_with(prefix_control)) / reads_total * 100 #divide number of control reads for each ASV by the number of total run reads (>> yields a percentage)
    control_props1 <- data.frame(row.names = row.names(control_props), #make empty dataframe to store values
                                 "prop" = c(rep(NA, nrow(control_props))))
    for(i in 1:nrow(control_props)){ #for each control:
      control_props1$prop[i] <- max(control_props[i, 1]:control_props[i, ncol(control_props)]) #take the maximum value of previously calculated percentage (in case there are several controls)
    }
    asvs_remove <- c(row.names(subset(control_props1, prop > 1))) #store ASV as removable if it occurs in the control at >1% of total run reads
    metbarc_clean1 <- metbarc_clean[!(row.names(metbarc_clean) %in% asvs_remove),]
  } else {
    metbarc_clean1 <- metbarc_clean
  }
  
  
  # Transpose dataframe for the next steps --------------------------
  ## Ditching non-samples columns and the columns representing control samples
  metbarc_clean12 <- metbarc_clean1
  metbarc_clean2 <- dplyr::select(metbarc_clean12, starts_with(prefix)) #selects sample columns
  dim(metbarc_clean2)
  
  ## Transform the bipartite matrices into a list (data-frame with links between predator and prey)
  metbarc.t <- data.frame(t(metbarc_clean2)) #transpose so I have individuals as ID variable (rownames)
  rownames(metbarc.t)
  colnames(metbarc.t)
  metbarc.t2 <- cbind.data.frame(reference=row.names(metbarc.t), metbarc.t) #creating a reference to use in melt function
  links <- melt(metbarc.t2, id.vars="reference",na.rm = T) #changes organization of data - uses individuals as ID
  head(links)
  colnames(links) <- c("predator", "prey", "weight")
  links[,1] = as.character(paste(links[,1]))
  links[,2] = as.character(paste(links[,2]))
  head(links)
  dim(links)
  
  ## Exclude all links with value 0 (zero) - need always to be removed due to melt function
  links1 <- subset(links, weight > 0)
  head(links1)
  dim(links1)
  
  ## Add column about number of reads of each ASV in the control
  links2 <- dplyr::select(metbarc_clean12, starts_with(prefix_control)) %>%  #selects control columns
    mutate(control_asvs = rowSums(.),
           prey = row.names(.)) %>% 
    dplyr::select(prey, control_asvs) %>% 
    left_join(links1, ., by = "prey") %>% as.data.frame
  attr(links2$control_asvs, "ATT") <- NULL #delete attribute that was inserted due to rowSums()
  
  ### Merging ASVs info 
  ## Subsetting original metabarcoding data to get taxa info
  asvs <- data.frame(asv = rownames(metbarc), subset(metbarc, select=c(phylum, class, order, family, genus, species)))
  head(asvs)
  
  ## Merging order with links
  links_asv <- merge(links2, asvs, by.x ="prey", by.y="asv", all.x = TRUE)
  head(links_asv)
  dim(links_asv)
  links_order <- links_asv %>%
    dplyr::select(c(predator, prey, weight, control_asvs, phylum, class, order, family, genus, species)) %>% #removes all other columns
    dplyr::rename(asv_phylum = phylum, asv_class = class, asv_order = order, 
                  asv_family = family, asv_genus = genus, asv_species = species)
  head(links_order)
  dim(links_order)
  
  ### Merging metabarcoding data with species and habitat data
  samples <- sample_list #importing database with faeces samples ID
  info <- subset(samples, select = c(lab.nbr, animal, species, site, season, year, landscape, collection_date, long, lat, sex, age))
  head(info)
  info$lab.nbr<-paste(prefix,info$lab.nbr,sep="") #add prefix to lab number to match sample database
  links_habitat <- merge(links_order, info, by.x ="predator", by.y="lab.nbr", all.x = TRUE)
  head(links_habitat)
  dim(links_habitat)
  
  
  # KEEP_CLASS --------------------------
  ## Only want to keep relevant prey taxa - only want to keep Arthropoda:Arachnida and Insecta for now
  levels(as.factor(links_habitat$asv_phylum))
  levels(as.factor(links_habitat$asv_class))
  if(is.null(keep_class)){
    links_clean <- links_habitat
  } else{
    links_clean <-links_habitat %>%
      filter(asv_class %in% keep_class) #default c("Arachnida","Insecta") - this can be changed to include other classes
  }
  dim(links_clean)
  head(links_clean)
  
  
  # HITS_CLEAN --------------------------
  links_hit <- data.frame()
  for (i in 1:length(unique(links_clean$predator))){
    hd <- links_clean[links_clean$predator == unique(links_clean$predator)[i],] #select links from a specific individual
    total_reads <- sum(hd$weight)
    hd <- hd %>% group_by(asv_phylum, asv_class, asv_order, asv_family, asv_genus, asv_species) %>% 
      summarize(reads_hit = sum(weight),
                reads_hit_proportion = reads_hit/total_reads)
    hd_join <- left_join(subset(links_clean, predator == unique(links_clean$predator)[i]), hd, 
                         by=c('asv_phylum', 'asv_class', 'asv_order', 'asv_family', 'asv_genus', 'asv_species'))
    links_hit <- rbind(links_hit, hd_join)
  } 
  head(links_hit)
  dim(links_hit)
  
  if(hits_clean == 0){
    links_clean1 <- links_hit
  }else if(hits_clean > 5){
    stop("Threshold too high (only values from 0 to 5 are accepted in asvs_clean)", call. = FALSE) 
  } else {
    check_deleted_hits <- links_hit %>% filter(reads_hit_proportion < (hits_clean/100)) %>% 
      group_by(asv_phylum, asv_class, asv_order, asv_family, asv_genus, asv_species) %>% 
      summarize(n_predators = n_distinct(predator),
                min_prop = min(reads_hit_proportion),
                max_prop = max(reads_hit_proportion)) %>% 
      arrange(asv_phylum, asv_class, asv_order, asv_family, asv_genus, asv_species)
    check_kept_hits <- links_hit %>% filter(reads_hit_proportion >= (hits_clean/100)) %>% 
      group_by(asv_phylum, asv_class, asv_order, asv_family, asv_genus, asv_species) %>% 
      summarize(n_predators = n_distinct(predator),
                min_prop = min(reads_hit_proportion),
                max_prop = max(reads_hit_proportion)) %>% 
      arrange(asv_phylum, asv_class, asv_order, asv_family, asv_genus, asv_species)
    
    ## Filter
    links_clean1 <- subset(links_hit, reads_hit_proportion >= (hits_clean/100))
  }
  dim(links_clean1)
  head(links_clean1)
  
  
  # ASVS_CLEAN --------------------------
  ### Excluding ASVs with less than 1% of total number of reads per sample
  ## Creating column with proportion of ASVs per individual
  links_propor <- data.frame()
  for (i in 1:length(unique(links_clean1$predator))){
    id <- links_clean1[links_clean1$predator == unique(links_clean1$predator)[i],] %>% #select links from a specific individual
      mutate(total_reads = sum(weight), #calculates total number of reads for that individual
             total_asvs = n_distinct(prey), #calculates total number of ASVs for that individual
             proportion = weight/total_reads) #calculates proportions of each ASV's reads to total reads for that individual
    links_propor <- rbind(links_propor, id)
  } 
  head(links_propor)
  dim(links_propor)

  if(asvs_clean == 0){
    links_clean2 <- links_propor
  }else if(asvs_clean > 5){
    stop("Threshold too high (only values from 0 to 5 are accepted in asvs_clean)", call. = FALSE) 
  } else {
    ## Filter
    links_clean2 <- subset(links_propor, proportion >= (asvs_clean/100))
  }
  dim(links_clean2)
  head(links_clean2)
  
  
  # CONTROL_ASVS_CLEAN --------------------------
  ### Remove ASVs from samples when they appear in controls at >1% (of the total occurrence of this ASV in the dataset.
  if(control_ASVs_clean == 0){
    links_clean3 <- links_clean2
  } else if(control_ASVs_clean > 0 & control_ASVs_clean < 100){
    links_clean3 <- links_clean2 %>% 
      
      ## For each ASV, add the number of reads in the control to the number of reads present in all samples/predators
      group_by(prey) %>% summarize(total_asv_reads = 
                                     sum(weight) + #sum of reads of a ASV in all samples
                                     unique(control_asvs)) %>% #unique, because the number of asv reads in the control is repeated for each sample/predator
      merge(links_clean2, .) %>% #add the new info to dataset
      ## Calculate per ASV the proportion of control reads to the total reads
      mutate(perc_control_of_total_asv_reads = control_asvs / total_asv_reads * 100,
             perc_control_of_total_asv_reads = replace_na(perc_control_of_total_asv_reads, 0)) #because NaN were created due to division with 0
    
    ## Checking which ASVs will be removed >> Check which cut-off makes sense for your dataset
    check_kept_df <- links_clean3 %>% 
      dplyr::select(total_asvs, total_asv_reads, control_asvs, perc_control_of_total_asv_reads, asv_phylum:asv_species) %>%  #make the dataset easier to read (slim down)
      filter(perc_control_of_total_asv_reads >= control_ASVs_clean)
    check_kept_df
    
    ## Remove ASVs that are found in the control at >1% (default) of their total occurrence in the dataset
    links_clean3 <- filter(links_clean3, perc_control_of_total_asv_reads < control_ASVs_clean) %>% 
      dplyr::select(-c(perc_control_of_total_asv_reads, total_asv_reads, control_asvs))
  } else {
    stop("Threshold too high (only values from 0 to 100 are accepted in control_ASVs_clean)", call. = FALSE) 
  }
  
  
  # REMOVE_NAORDERS --------------------------
  ### Checking order to decide if I should pull only the records that got a hit to ORDER
  if(remove_NAorders == T){
    sum(is.na(links_clean3$asv_order))
    levels(as.factor(links_clean3$asv_order))
    links_clean4 <- subset(links_clean3, asv_order != "NA")
  } else {
    links_clean4 <- links_clean3
  }
  dim(links_clean4)
  
  
  # REMOVE_NAFAMILY --------------------------
  ### Checking family to decide if I should pull only the records that got a hit to family
  if(remove_NAfamily == T){
    sum(is.na(links_clean4$asv_family))
    levels(as.factor(links_clean4$asv_family))
    links_clean5 <- subset(links_clean4, asv_family != "NA")
  } else {
    links_clean5 <- links_clean4
  }
  dim(links_clean5)
  
  
  # DESIRED_SPECIES --------------------------
  ### Only keeping the predator species that I want 
  if(is.null(desired_species)){
    links_filter <- links_clean5
  } else{
    links_filter <- links_clean5 %>%
      filter(species %in% desired_species) #e.g. c("Hipposideros ruber")
  }
  head(links_filter)
  dim(links_filter)
  
  
  # update total_reads, total_asvs, proportions --------------------------
  links_filter1 <- links_filter %>% 
    dplyr::select(-c(total_reads,total_asvs,proportion,reads_hit,reads_hit_proportion))
  
  links_filter2 <- group_by(links_filter1, predator) %>% summarize(total_reads = sum(weight),
                                                                   total_asvs = length(unique(prey))) %>% 
    left_join(links_filter1, ., by = "predator") %>% 
    mutate(proportion = weight / total_reads)
  links_filter3 <- links_filter2 %>%
    group_by(predator, asv_phylum, asv_class, asv_order, asv_family, asv_genus, asv_species) %>% reframe(reads_hit = sum(weight),
                                                                                                 reads_hit_proportion = reads_hit/total_reads) %>% unique %>% 
    left_join(links_filter2, ., by = c("predator", "asv_phylum", "asv_class", "asv_order", "asv_family", "asv_genus", "asv_species")) %>%
    dplyr::select(where(~ !(all(is.na(.)) | all(. == "")))) #remove all empty columns that were added to make the code work
  

  # Final dataset --------------------------
  return(links_filter3)
}

