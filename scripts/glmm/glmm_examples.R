#' ---
#' title: Examples of GLMMs with diet metabarcoding data
#' author: Crinan Jarrett
#' date: 24 Sept 2021
#' ---

library(lme4)
library(reshape2)

zbj <- read.csv("data/metabarcoding/zbj_consensus.tsv", sep = "\t",row.names = 1,header = TRUE, na.strings=c("NA", "NULL", "", ".")) 
colnames(zbj)

###Reading my function that returns final data frame with all data organized and cleaned - see generate_final_data.r for more details
source('scripts/generate_final_data.r')

zbj_data <- final_metbar(data = zbj,remove_samples=T,otus_clean=T,keep_class=c("Arachnida","Insecta"),remove_NAorders=T,remove_NAfamily=F,desired_species=NULL)

###Casting dataframe into matrix predator x prey
data_recast<-dcast(zbj_data, formula=predator+total+landscape~prey,value.var = "weight")

###Converting NAs into 0 for presence absence
data_occurrence<-data_recast
data_occurrence[,4:50]<-sapply(data_occurrence[,4:50],as.numeric)
data_occurrence <- data_occurrence %>% dplyr::mutate_if(is.numeric,~replace(., is.na(.), 0))
data_occurrence <- data_occurrence %>% dplyr::mutate_if(is.numeric, ~replace(., . > 0, 1)) 

## Scenario: faecal samples from bats in 4 different lanscapes (on a gradient of degradation)
## Explanatory variable: Landscape (categorical)

## FOO

## %FOO

## Weighted Occurrence

## RRA

## Number of reads