
##################################################
## Project: How to metabarcode (Rachel et al 2024)
## Script purpose: Examples of GLMs with diet metabarcoding data 
## Date: 29.03.2024
## Author: Crinan Jarrett (crinan.jarrett@gmail.com )
## Packages: lme4, reshape2, DHARMa, MASS and openxlsx 
## Notes: you will need the function in organize&clean_metabarcoding.r to run this script 
###################################################

library(lme4)
library(reshape2)
library(DHARMa)
library(MASS)
library(openxlsx)
library(glmmTMB)

## Scenario: faecal samples from bats in 4 different landscapes and on a gradient of management
## Explanatory variables: Management (continuous) and Landscape (categorical)

### 1. FORMATTING DATA INTO DATAFRAME SUITABLE FOR GLMs
# importing metabarcoding data 
zbj <- read.csv2("./data/16Farms_ASV_Table.csv", dec = ".", row.names = 1,header = TRUE, check.names=FALSE,na.strings=c("NA", "NULL", "", ".")) 

# importing species information for each sample
samples_list <- read.xlsx("./data/faeces_sample_database.xlsx", sheet = "samples")

# Reading function that returns final data frame with all data organized and cleaned - see organize&clean_metabarcoding.r for more details
source('./scripts/1.organize&clean_metabarcoding.R')
zbj_data <- final_metbar(data = zbj,sample_list = samples_list,prefix_control = "Control",
                         remove_samples=F,remove_control_asvs=1,asvs_clean=1,
                         control_asvs_clean=1, hits_clean=0,keep_class=c("Arachnida","Insecta"),
                         remove_NAorders=T,remove_NAfamily=F,desired_species=NULL)

zbj_data$site<-as.factor(zbj_data$site)

# Simulate covariate data (for demonstration of multiple variables)
management<-data.frame(site=levels(zbj_data$site),
                       manage=rnorm(length(levels(zbj_data$site)),0,1))

# now merge covariate data with diet dataframe
zbj_data<-merge(zbj_data,management,by="site",
                all.x=TRUE)

# group dataframe by ASV order as we will use this taxonomic level for analyses
zbj_data_order<- zbj_data %>%
  group_by(site,predator,ASV_order,year,season,landscape,species,total_reads,manage) %>%
  summarise(order_asvs=n(),sum_weight=sum(weight))

# calculate total ASVs per sample (used later for Weighted Abundance)
zbj_total<- zbj_data %>%
  group_by(predator) %>%
  summarise(total_asvs=n())

# merge with order-level grouped dataframe
zbj_data_order<-merge(zbj_data_order,zbj_total,by="predator")

# casting dataframe into matrix predator x prey
data_recast<-dcast(zbj_data_order, formula=predator+manage+total_reads+landscape~ASV_order,value.var = "sum_weight")

## For occurrence models
# converting NAs into 0  and >1 into 1 for presence absence
data_occurrence<-data_recast 
data_occurrence[,5:ncol(data_occurrence)] <- data_occurrence[,5:ncol(data_occurrence)] %>% dplyr::mutate(replace(., is.na(.), 0))
data_occurrence[,5:ncol(data_occurrence)] <- data_occurrence[,5:ncol(data_occurrence)] %>% dplyr::mutate(replace(., . > 0, 1)) 

# make landscape into factor so it can be used as grouping variable
data_occurrence$landscape<-as.factor(data_occurrence$landscape)

# group occurrence data by landscape and management (explanatory variables)
occurrence_grouped<- data_occurrence %>%
  group_by(landscape,manage) %>%
  dplyr::summarise(sample_n=n(),across("Coleoptera":"Lepidoptera",sum))

### 2. GENERALISED LINEAR MODELS

## We will use Diptera as the target group for this analysis as it shows more variation than Lepidoptera across the landscapes.

## EXAMPLE 1: FOO (Frequency of occurrence; number of individuals with Diptera present in diet)
mod_foo<-glm(Diptera ~ manage + landscape, family = "poisson", data = occurrence_grouped)
summary(mod_foo) #the summary shows that landscape has an effect on FOO of Diptera, with landscape 'Bokito' having highest FOO, followed by 'Ayos' (baseline) and lowest in 'Konye'  
DHARMa::testDispersion(mod_foo,alternative="greater") #data not overdispersed
DHARMa::testDispersion(mod_foo,alternative="less") #data not underdispersed

plot(occurrence_grouped$Diptera~occurrence_grouped$landscape)
plot(occurrence_grouped$Diptera~occurrence_grouped$manage)

## EXAMPLE 2: %FOO (% Frequency of occurrence; proportion of individuals with Diptera present in diet)
mod_percfoo<-glm(cbind(Diptera,sample_n) ~ manage + landscape, family = "binomial", data = occurrence_grouped)
summary(mod_percfoo) #here the summary indicates that management doesn't have an effect on %FOO but landscape does 
DHARMa::testDispersion(mod_percfoo,alternative="greater") #data not overdispersed
DHARMa::testDispersion(mod_percfoo,alternative="less") #data not underdispersed

plot((occurrence_grouped$Diptera/occurrence_grouped$sample_n)~occurrence_grouped$landscape)
plot((occurrence_grouped$Diptera/occurrence_grouped$sample_n)~occurrence_grouped$manage)

## EXAMPLE 3: Weighted Occurrence (number of detected ASVs corresponding to Diptera divided by the total number of ASVs detected in each individual)
# first, we need to filter dataframe to select only Diptera rows
Diptera_data_order<- zbj_data_order %>%
  filter(ASV_order=="Diptera")
Diptera_data_order$landscape<-as.factor(Diptera_data_order$landscape)

mod_wocc<-glm(cbind(order_asvs,total_asvs) ~ manage + landscape, family = "binomial", data = Diptera_data_order)
summary(mod_wocc) #here the summary indicates that neither landscape or management significantly affect WOcc of Diptera
DHARMa::testDispersion(mod_wocc,alternative="greater") #data not overdispersed
DHARMa::testDispersion(mod_wocc,alternative="less") #data underdispersed, so we can try alternative distribution family. See below for potential solution: log-transforming values and then using Gaussian distribution with offset

mod_wocc_log<-glm(log(order_asvs) ~ manage + landscape, offset = log(total_asvs), data = Diptera_data_order)
summary(mod_wocc_log) #indicates that the WOcc is not significantly different. 
DHARMa::testDispersion(mod_wocc_log,alternative="greater") #data not overdispersed
DHARMa::testDispersion(mod_wocc_log,alternative="less") #The not underdispersed

plot((Diptera_data_order$order_asvs/Diptera_data_order$total_asvs)~Diptera_data_order$landscape)
plot((Diptera_data_order$order_asvs/Diptera_data_order$total_asvs)~Diptera_data_order$manage)

## EXAMPLE 4: RRA (Relative Read Abundance; proportion of reads represented by Diptera)
mod_rra<-glm(cbind(sum_weight,total_reads) ~ manage + landscape, family = "binomial", data = Diptera_data_order)
summary(mod_rra) #summary seems to indicate that both management and landscape have significant effect on RRA of Diptera, however, check dispersion (below)
DHARMa::testDispersion(mod_rra,alternative="greater") #data overdispersed, suggest trying alternative distribution family. Below we test log-transforming values and then using Gaussian distribution with offset
DHARMa::testDispersion(mod_rra,alternative="less") #data not underdispersed

mod_rra.log<-glm(log(sum_weight) ~ manage + landscape, offset = log(total_reads), data = Diptera_data_order)
summary(mod_rra.log) #model output now suggests that neither management nor landscape have significant effects on RRA (this is more what we would expect based on plots below)
DHARMa::testDispersion(mod_rra.log,alternative="greater") #data no longer overdispersed
DHARMa::testDispersion(mod_rra.log,alternative="less") #data not underdispersed

plot((Diptera_data_order$sum_weight/Diptera_data_order$total_reads)~Diptera_data_order$landscape)
plot((Diptera_data_order$sum_weight/Diptera_data_order$total_reads)~Diptera_data_order$manage)

## EXAMPLE 5: Number of reads
mod_reads<-glm(sum_weight ~ manage + landscape, family = "poisson", data = Diptera_data_order)
summary(mod_reads) #summary seems to indicate that both management and landscape have significant effect on number of reads of Diptera, however, check dispersion (below)
DHARMa::testDispersion(mod_reads,alternative="greater") #data overdispersed, suggest trying alternative distribution family. Below we test negative binomial
DHARMa::testDispersion(mod_reads,alternative="less") #data not underdispersed

mod_reads.nb<-glm.nb(sum_weight ~ manage + landscape, data = Diptera_data_order)
summary(mod_reads.nb) #model output still suggests that landscape and (marginally) management have an effect on number of reads. 
DHARMa::testDispersion(mod_reads.nb,alternative="greater") #The overdispersion has been improved, though still could cause issues
DHARMa::testDispersion(mod_reads.nb,alternative="less") #data not underdispersed

plot(Diptera_data_order$sum_weight~Diptera_data_order$landscape)
plot(Diptera_data_order$sum_weight~Diptera_data_order$manage)

### 3. MULTIVARIATE GENERALISED LINEAR MODELS: model joint distribution of prey groups as a function of explanatory variable (landscape)
library(mvabund)

# select columns that correspond to 'abundance' (i.e., number of reads) of each prey group in each sample
reads_multivar<- data_recast[,5:ncol(data_recast)]
# replace NAs (absences) with 0s
reads_multivar[is.na(reads_multivar)]<-0
# generate mvabund object, used in functions from mvabund package
reads_dat<-mvabund(reads_multivar)
# generate vector of explanatory variable (landscape)
landscape<-data_recast[,4]

# fit multivariate model. We selected negative binomial family to account for overdispersion in read abundance data.
mod_multivar<-manyglm(reads_dat~landscape,family="negative.binomial")
mod_multivar
summary(mod_multivar) #output shows a that there is a marginally significant effect of landscape on the composition of prey groups in diet

