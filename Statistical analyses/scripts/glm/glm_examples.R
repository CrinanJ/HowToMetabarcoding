#' ---
#' title: Examples of GLMs with diet metabarcoding data
#' author: Crinan Jarrett
#' date: 21 Jan 2022
#' ---
 
library(lme4)
library(reshape2)
library(DHARMa)
library(MASS)

## Scenario: faecal samples from bats in 4 different landscapes and on a gradient of management
## Explanatory variables: Management (continuous) and Landscape (categorical)

### 1. FORMATTING DATA INTO DATAFRAME SUITABLE FOR GLMs
# importing metabarcoding data 
zbj <- read.csv2("data/16Farms_OTU_Table_vsearch.csv", dec = ".", row.names = 1,header = TRUE, check.names=FALSE,na.strings=c("NA", "NULL", "", ".")) 

# importing species information for each sample
library(openxlsx)#read excel and sheet
samples_list <- read.xlsx("data/faeces_sample_database.xlsx", sheet = "samples")

# Reading function that returns final data frame with all data organized and cleaned - see organize&clean_metabarcoding.r for more details
source('scripts/organize&clean_metabarcoding.r')
zbj_data <- final_metbar(data = zbj,sample_list = samples_list, remove_samples=F,otus_clean=1, keep_class=c("Arachnida","Insecta"),remove_NAorders=T,remove_NAfamily=F,desired_species=NULL)
zbj_data$farm<-as.factor(zbj_data$farm)

# Simulate covariate data
management<-data.frame(farm=levels(zbj_data$farm),
                       manage=rnorm(length(levels(zbj_data$farm)),0,1))

# now merge covariate data with diet dataframe
zbj_data<-merge(zbj_data,management,by="farm",
                all.x=TRUE)

# group dataframe by OTU order as we will use this taxonomic level for analyses
zbj_data_order<- zbj_data %>%
  group_by(farm,predator,otu_order,year,season,landscape,species,total_reads,manage) %>%
  summarise(order_otus=n(),sum_weight=sum(weight))

# calculate total OTUs per sample (used later for Weighted Abundance)
zbj_total<- zbj_data %>%
  group_by(predator) %>%
  summarise(total_otus=n())

# merge with order-level grouped dataframe
zbj_data_order<-merge(zbj_data_order,zbj_total,by="predator")

# casting dataframe into matrix predator x prey
data_recast<-dcast(zbj_data_order, formula=predator+manage+total_reads+landscape~otu_order,value.var = "sum_weight")

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

## We will use Lepidoptera are target group for this analysis

## EXAMPLE 1: FOO (Frequency of occurrence; number of individuals with Lepidoptera present in diet)
mod_foo<-glm(Lepidoptera ~ manage + landscape, family = "poisson", data = occurrence_grouped)
summary(mod_foo) #the summary shows that landscape has an effect on FOO of Lepidoptera, with landscape 'Bokito' having highest FOO, followed by 'Ayos' (baseline) and finally 'Konye'  
DHARMa::testDispersion(mod_foo,alternative="greater") #data not overdispersed
DHARMa::testDispersion(mod_foo,alternative="less") #data not underdispersed

plot(occurrence_grouped$Lepidoptera~occurrence_grouped$landscape)
plot(occurrence_grouped$Lepidoptera~occurrence_grouped$manage)

## EXAMPLE 2: %FOO (% Frequency of occurrence; proportion of individuals with Lepidoptera present in diet)
mod_percfoo<-glm(cbind(Lepidoptera,sample_n) ~ manage + landscape, family = "binomial", data = occurrence_grouped)
summary(mod_percfoo) #here the summary indicates that neither landscape or management significantly affect %FOO of Lepidoptera
DHARMa::testDispersion(mod_percfoo,alternative="greater") #data not overdispersed
DHARMa::testDispersion(mod_percfoo,alternative="less") #data not underdispersed

plot((occurrence_grouped$Lepidoptera/occurrence_grouped$sample_n)~occurrence_grouped$landscape)
plot((occurrence_grouped$Lepidoptera/occurrence_grouped$sample_n)~occurrence_grouped$manage)

## EXAMPLE 3: Weighted Occurrence (number of detected MOTUs corresponding to Lepidoptera divided by the total number of MOTUs detected in each individual)
# first, we need to filter dataframe to select only Lepidoptera rows
lep_data_order<- zbj_data_order %>%
  filter(otu_order=="Lepidoptera")
lep_data_order$landscape<-as.factor(lep_data_order$landscape)

mod_wocc<-glm(cbind(order_otus,total_otus) ~ manage + landscape, family = "binomial", data = lep_data_order)
summary(mod_wocc) #here the summary indicates that neither landscape or management significantly affect WOcc of Lepidoptera
DHARMa::testDispersion(mod_wocc,alternative="greater") #data not overdispersed
DHARMa::testDispersion(mod_wocc,alternative="less") #data not underdispersed

plot((lep_data_order$order_otus/lep_data_order$total_otus)~lep_data_order$landscape)
plot((lep_data_order$order_otus/lep_data_order$total_otus)~lep_data_order$manage)

## EXAMPLE 4: RRA (Relative Read Abundance; proportion of reads represented by Lepidoptera)
mod_rra<-glm(cbind(sum_weight,total_reads) ~ manage + landscape, family = "binomial", data = lep_data_order)
summary(mod_rra) #summary seems to indicate that both management and landscape have significant effect on RRA of Lepidoptera, however, check dispersion (below)
DHARMa::testDispersion(mod_rra,alternative="greater") #data overdispersed, suggest trying alternative distribution family. Below we test log-transforming values and then using Gaussian distribution with offset
DHARMa::testDispersion(mod_rra,alternative="less") #data not underdispersed

mod_rra.log<-glm(log(sum_weight) ~ manage + landscape, offset = log(total_reads), data = lep_data_order)
summary(mod_rra.log) #model output now suggests no significant effects of management or landscape on RRA (this is more what we would expect based on plots below)
DHARMa::testDispersion(mod_rra.log,alternative="greater") #data no longer overdispersed
DHARMa::testDispersion(mod_rra.log,alternative="less") #data not underdispersed

plot((lep_data_order$sum_weight/lep_data_order$total_reads)~lep_data_order$landscape)
plot((lep_data_order$sum_weight/lep_data_order$total_reads)~lep_data_order$manage)

## EXAMPLE 5: Number of reads
mod_reads<-glm(sum_weight ~ manage + landscape, family = "poisson", data = lep_data_order)
summary(mod_reads) #summary seems to indicate that both management and landscape have significant effect on number of reads of Lepidoptera, however, check dispersion (below)
DHARMa::testDispersion(mod_reads,alternative="greater") #data overdispersed, suggest trying alternative distribution family. Below we test negative binomial
DHARMa::testDispersion(mod_reads,alternative="less") #data not underdispersed

mod_reads.nb<-glm.nb(sum_weight ~ manage + landscape, data = lep_data_order)
summary(mod_reads.nb) #model output now suggests a significant effect of management (but not landscape) on number of reads. This is what we would expect based on plots below
DHARMa::testDispersion(mod_reads.nb,alternative="greater") #data still overdispersed but less so
DHARMa::testDispersion(mod_reads.nb,alternative="less") #data not underdispersed

plot(lep_data_order$sum_weight~lep_data_order$landscape)
plot(lep_data_order$sum_weight~lep_data_order$manage)
