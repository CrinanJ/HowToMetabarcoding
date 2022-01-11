#' ---
#' title: Examples of GLMMs with diet metabarcoding data
#' author: Crinan Jarrett
#' date: 24 Sept 2021
#' ---

library(lme4)
library(reshape2)

###importing metabarcoding data 
zbj <- read.csv2("data/16Farms_OTU_Table_withConfidence.csv", dec = ".", row.names = 1,header = TRUE, check.names=FALSE,na.strings=c("NA", "NULL", "", ".")) 
colnames(zbj)

###importing species information for each sample
library(openxlsx)#read excel and sheet
samples_list <- read.xlsx("data/feces_sample_database.xlsx", sheet = "samples")
colnames(samples_list)

###Reading function that returns final data frame with all data organized and cleaned - see organize&clean_metabarcoding.r for more details
source('scripts/organize&clean_metabarcoding.r')

zbj_data <- final_metbar(data = zbj,sample_list = samples_list, remove_samples=F,otus_clean=T,keep_class=NULL,remove_NAorders=T,remove_NAfamily=F,desired_species=NULL)
colnames(zbj_data)

###Casting dataframe into matrix predator x prey
data_recast<-dcast(zbj_data, formula=predator+total_reads+landscape~prey,value.var = "weight")

###Converting NAs into 0 for presence absence
data_occurrence<-data_recast
data_occurrence[,4:50]<-sapply(data_occurrence[,4:50],as.numeric)
data_occurrence <- data_occurrence %>% dplyr::mutate_if(is.numeric,~replace(., is.na(.), 0))
data_occurrence <- data_occurrence %>% dplyr::mutate_if(is.numeric, ~replace(., . > 0, 1)) 

## Scenario: faecal samples from bats in 4 different landscapes (on a gradient of degradation)
## Explanatory variable: Landscape (categorical)

## FOO

## %FOO

## Weighted Occurrence

## RRA

## Number of reads