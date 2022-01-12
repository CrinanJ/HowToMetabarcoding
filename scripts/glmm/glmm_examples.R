#' ---
#' title: Examples of GLMMs with diet metabarcoding data
#' author: Crinan Jarrett
#' date: 24 Sept 2021
#' ---
 
library(lme4)
library(reshape2)
library(DHARMa)

###importing metabarcoding data 
zbj <- read.csv2("data/16Farms_OTU_Table_withConfidence.csv", dec = ".", row.names = 1,header = TRUE, check.names=FALSE,na.strings=c("NA", "NULL", "", ".")) 
colnames(zbj)
rownames(zbj)<-paste("otu",seq(1,nrow(zbj)),sep="_")

###importing species information for each sample
library(openxlsx)#read excel and sheet
samples_list <- read.xlsx("data/faeces_sample_database.xlsx", sheet = "samples")
colnames(samples_list)

###Reading function that returns final data frame with all data organized and cleaned - see organize&clean_metabarcoding.r for more details
source('scripts/organize&clean_metabarcoding.r')
zbj_data <- final_metbar(data = zbj,sample_list = samples_list, remove_samples=F,otus_clean=T,keep_class=NULL,remove_NAorders=T,remove_NAfamily=F,desired_species=NULL)
colnames(zbj_data)

###Casting dataframe into matrix predator x prey
data_recast<-dcast(zbj_data, formula=predator+total_reads+landscape~prey,value.var = "weight")

###Converting NAs into 0 for presence absence
data_occurrence<-data_recast 
#data_occurrence[,4:ncol(data_occurrence)]<-sapply(data_occurrence[,ncol(data_occurrence)],as.numeric)
data_occurrence[,4:ncol(data_occurrence)] <- data_occurrence[,4:ncol(data_occurrence)] %>% dplyr::mutate(replace(., is.na(.), 0))
data_occurrence[,4:ncol(data_occurrence)] <- data_occurrence[,4:ncol(data_occurrence)] %>% dplyr::mutate(replace(., . > 0, 1)) 
data_occurrence$landscape<-as.factor(data_occurrence$landscape)

# group occurrence data by landscape
occurrence_grouped<- data_occurrence %>%
  group_by(landscape) %>%
  dplyr::summarise(across("otu_1":"otu_998",sum))

## Scenario: faecal samples from bats in 4 different landscapes (on a gradient of degradation)
## Explanatory variable: Landscape (categorical)

## FOO: Frequency of occurrence (number of individuals with certain OTU present in diet)
mod_foo<-glm(otu_6 ~ landscape, family = "poisson", data = occurrence_grouped)
summary(mod_foo)
DHARMa::testOverdispersion(mod_foo)
  
## %FOO

## Weighted Occurrence

## RRA

## Number of reads