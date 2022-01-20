#' ---
#' title: Examples of GLMs with diet metabarcoding data
#' author: Crinan Jarrett
#' date: 13 Jan 2022
#' ---
 
library(lme4)
library(reshape2)
library(DHARMa)
library(MASS)

## Scenario: faecal samples from bats in 4 different landscapes and on a gradient of management
## Explanatory variables: Management (continuous) and Landscape (categorical)

### 1. FORMATTING DATA INTO DATAFRAME SUITABLE FOR GLMs
# importing metabarcoding data 
zbj <- read.csv2("data/16Farms_OTU_Table_withConfidence.csv", dec = ".", row.names = 1,header = TRUE, check.names=FALSE,na.strings=c("NA", "NULL", "", ".")) 

# making rownames easier to handle
rownames(zbj)<-paste("otu",seq(1,nrow(zbj)),sep="_")

# importing species information for each sample
library(openxlsx)#read excel and sheet
samples_list <- read.xlsx("data/faeces_sample_database.xlsx", sheet = "samples")

# Reading function that returns final data frame with all data organized and cleaned - see organize&clean_metabarcoding.r for more details
source('scripts/organize&clean_metabarcoding.r')
zbj_data <- final_metbar(data = zbj,sample_list = samples_list, remove_samples=F,otus_clean=1,keep_class=NULL,remove_NAorders=T,remove_NAfamily=F,desired_species=NULL)
zbj_data$farm<-as.factor(zbj_data$farm)

# Simulate covariate data
management<-data.frame(farm=levels(zbj_data$farm),
                       manage=rnorm(17,0,1))

# now merge covariate data with diet dataframe
zbj_data<-merge(zbj_data,management,by="farm",
                all.x=TRUE)

# group dataframe by OTU order as we will you this taxonomic level for analyses
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

## EXAMPLE 1: FOO (Frequency of occurrence; number of individuals with certain OTU present in diet)
mod_foo<-glm(Coleoptera ~ manage + landscape, family = "poisson", data = occurrence_grouped)
summary(mod_foo) #the summary shows that landscape has an effect on FOO of Lepidotera, with landscape 'Bokito' having highest FOO, followed by 'Ayos' (baseline) and finally 'Konye'  
DHARMa::testOverdispersion(mod_foo) #not overdispersed so Poisson is ok

plot(occurrence_grouped$Coleoptera~occurrence_grouped$landscape)
plot(occurrence_grouped$Coleoptera~occurrence_grouped$manage)

## EXAMPLE 2: %FOO (% Frequency of occurrence; proportion of individuals with certain OTU present in diet)
mod_percfoo<-glm(cbind(Coleoptera,sample_n) ~ manage + landscape, family = "binomial", data = occurrence_grouped)
summary(mod_percfoo) #here the summary indicates that neither landscape or management significantly affect %FOO of OTU_6
DHARMa::testDispersion(mod_percfoo) # no overdispersion

plot((occurrence_grouped$Coleoptera/occurrence_grouped$sample_n)~occurrence_grouped$landscape)
plot((occurrence_grouped$Coleoptera/occurrence_grouped$sample_n)~occurrence_grouped$manage)

## EXAMPLE 3: Weighted Occurrence (number of detected MOTUs corresponding to the target prey group divided by the total number of MOTUs detected in each individual)
# first, we need to filter dataframe to select only Coleoptera rows
col_data_order<- zbj_data_order %>%
  filter(otu_order=="Coleoptera")
col_data_order$landscape<-as.factor(col_data_order$landscape)

mod_wocc<-glm(cbind(order_otus,total_otus) ~ manage + landscape, family = "binomial", data = col_data_order)
summary(mod_wocc)
DHARMa::testDispersion(mod_wocc) #test shows underdispersion?

plot((col_data_order$order_otus/col_data_order$total_otus)~col_data_order$landscape)
plot((col_data_order$order_otus/col_data_order$total_otus)~col_data_order$manage)

## EXAMPLE 4: RRA (Relative Read Abundance; proportion of reads represented by certain OTU or OTU group)
mod_rra<-glm(cbind(sum_weight,total_reads) ~ manage + landscape, family = "binomial", data = col_data_order)
summary(mod_rra)
DHARMa::testDispersion(mod_rra) # not overdispersed

plot((col_data_order$sum_weight/col_data_order$total_reads)~col_data_order$landscape)
plot((col_data_order$sum_weight/col_data_order$total_reads)~col_data_order$manage)

## EXAMPLE 5: Number of reads
mod_reads<-glm(sum_weight ~ manage + landscape, family = "poisson", data = col_data_order)
summary(mod_reads)
DHARMa::testDispersion(mod_reads) #data very overdispersed, so try negative binomial family

mod_reads.nb<-glm.nb(sum_weight ~ manage + landscape, data = col_data_order)
summary(mod_reads.nb)
DHARMa::testOverdispersion(mod_reads.nb) #no longer overdispersed

plot(col_data_order$sum_weight~col_data_order$landscape)
plot(col_data_order$sum_weight~col_data_order$manage)
