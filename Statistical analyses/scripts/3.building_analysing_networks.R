
##################################################
## Project: How to metabarcode (Rachel et al 2022)
## Script purpose: creating, plotting and analyzing bipartite networks 
## Date: 17/03/2022
## Author: Diogo F. Ferreira (ferreiradfa@gmail.com)
## Packages: openxlsx, iNEXT, bipartite, RColorBrewer, GISTools and ggplot2 
## Notes: you will need the function to generate_final_data.r to run this script
###################################################

###importing metabarcoding data 
#zbj primer: 224 samples, 2 controls and 2889 unique ASVS
zbj <- read.csv2("./data/16Farms_ASV_Table.csv", dec = ".", row.names = 1,header = TRUE, check.names=FALSE,na.strings=c("NA", "NULL", "", ".")) 
dim(zbj)
colnames(zbj)
head(zbj)

###importing species information for each sample
library(openxlsx)#read excel and sheet
samples_list <- read.xlsx("./data/faeces_sample_database.xlsx", sheet = "samples")
colnames(samples_list)
head(samples_list)

###Reading function that returns final data frame with all data organized and cleaned - see organise&clean_metabarcoding.r for more details
source('./scripts/1.organize&clean_metabarcoding.r')

unique(zbj$class)
#remove all ASVs that represent less than 1%, keeping only Arachnida and Insecta, and remove ASVs not identified until order
zbj_data <- final_metbar(data = zbj,sample_list = samples_list,prefix_control = "Control",
                         remove_samples=F,remove_control_asvs=1,asvs_clean=1,
                         control_asvs_clean=1, hits_clean=0,keep_class=c("Arachnida","Insecta"),
                         remove_NAorders=T,remove_NAfamily=F,desired_species=NULL)

head(zbj_data)
dim(zbj_data)
#2067  rows and 20 columns
n_distinct(zbj_data$prey)#number of different asvs
#1217
n_distinct(zbj_data$predator)#number of samples 
#265


####aggregate data by landscape and species
colnames(zbj_data)#keeping only species, prey, asv_order, landscape, weight, proportion and creating column for asv FOO
zbj_aggr <- zbj_data %>%
  group_by(predator=species,prey,ASV_order,landscape)%>%
  summarise(weight=sum(weight),proportion=mean(proportion),FOO=n())
head(zbj_aggr)
dim(zbj_aggr)
#1749 rows and 7 columns
levels(as.factor(zbj_aggr$predator))#check how and what bats/bird species I have
#11 

##checking how many samples per species and landscape
zbj_species <- zbj_data %>%
  group_by(landscape,animal,species)%>%
  summarise(samples = n_distinct(predator))%>%
  mutate(total_samples = sum(samples),n_species = n_distinct(species))%>%#creates columns with total samples per landscaoe and total number of species
  rename("predator"=species)
head(zbj_species)

#saving list of samples per species 
write.csv(zbj_species,"n_samples_network_zbj.csv",row.names = FALSE)

##I'll remove species that have less than 5 samples (ideally higher the cut-off the better) per landscape 
species_remove_zbj <- unique(zbj_species[zbj_species$samples < 5,])
levels(as.factor(species_remove_zbj$predator))


zbj_filtered <- data.frame()
for (i in 1:length(unique(zbj_aggr$landscape))){
  land <- zbj_aggr[zbj_aggr$landscape==unique(zbj_aggr$landscape)[i],]#select rows from a landscape
  remove <- species_remove_zbj[species_remove_zbj$landscape == unique(zbj_aggr$landscape)[i],]$predator
  final <- land[!land$predator %in% remove,] #removes species with only one sample
  zbj_filtered <- rbind(zbj_filtered, final)
}
head(zbj_filtered)
dim(zbj_filtered)
#1718 rows and 7 columns
levels(as.factor(zbj_filtered$predator))
#10 species

####standardizing links (edges) by number of samples
data_network <- merge(zbj_filtered, zbj_species, by = c("predator","landscape"))
head(data_network)

data_network$FOO_std <- data_network$FOO/data_network$samples*100
head(data_network)


#####################Checking sample coverage#####################
##creating a empty data frame to store INEXT values
#creating a matrix for network
library(reshape2)
#creating matrix for inext
matrix_inext <- dcast(data_network, prey ~ landscape, length)
head(matrix_inext)

#adding total row
matrix_inext <- rbind(data.frame("prey" = "Total", "Ayos"=sum(matrix_inext$Ayos),"Bokito"=sum(matrix_inext$Bokito),"Konye"=sum(matrix_inext$Konye)),matrix_inext)
head(matrix_inext)

#give names to rows
rownames(matrix_inext) <- matrix_inext$prey
head(matrix_inext)

#removing ID column
matrix_inext <- matrix_inext[-1]
head(matrix_inext)


##calculating rarefaction curves
library(iNEXT)
View(matrix_inext)
icurves <- iNEXT(matrix_inext, datatype="incidence_freq", nboot=100)
icurves

##Save rarefaction curves results##
write.csv(icurves$iNextEst,"./outputs/results/network analysis/icurves.estimates.csv")
write.csv(icurves$AsyEst, "./outputs/results/network analysis/icurves.asymptote.csv")

##Plot rarefaction curves##
igraph_diversity<-ggiNEXT(icurves, type = 1)#1 for species diversity, 2 for sample coverage and 3 for species diversity vs sample coverage (gives samples completeness) 
igraph_diversity #ASVs richness is still far from total sample coverage. Our network will not show a full representation of the diet. Also, (some) differences in coverage between landscapes. This will limit the comparsions that we can do between landscapes. However, we will used null networks to minimize the effect of this. 

igraph_coverage<-ggiNEXT(icurves, type = 2)
igraph_coverage# to have a good representation of the diet we would need a higher number of samples (e.g. 800 samples to reach only 0.5 of coverage)

##estimating richness at certain value of sampling 
estimateD(matrix_inext, datatype="incidence_freq", base="coverage", level=0.8) # to know how many samples we need to reach 80% coverage in each landscape



##############Bipartite network - on bipartite package##############
####generating an network based on a adjacency matrix and creating an each network based on landscape and using FOO_std)
#if you want to generate networks for other variables (e.g., season), you need to change webID
library(bipartite)
zbj_nets <- frame2webs(data.frame(lower=data_network$prey,higher=data_network$predator,webID=data_network$landscape,freq=data_network$FOO_std),type.out="list")
lapply(zbj_nets, dim)#checking dimensions of all data frames

unlist(lapply(zbj_nets, function(x) sum(x>=1)))#number of links per network
#Ayos Bokito  Konye 
#597  495     626 

unlist(lapply(zbj_nets, sum))#total asvs frequency per network
#Ayos         Bokito       Konye 
#5148.716    4308.164    5967.342  

unlist(lapply(zbj_nets, nrow))#number of asvs per network
#Ayos Bokito  Konye 
#525   382    513 

unlist(lapply(zbj_nets, ncol))#number of species per network
#Ayos Bokito  Konye 
#7      5      8 

#####This tiny chunk shows how many asvs are shared between networks
samples <- lapply(zbj_nets, rownames)
Reduce(intersect, samples)
length(Reduce(intersect, samples))
#37 asvs shared
(length(Reduce(intersect, samples))/length(unique(zbj_filtered$prey)))*100
#3.070539% of asvs are shared between all networks


####Plotting bipartite network
##creating pallet (28 colors and colorblind) to give to each arthropod order
library(RColorBrewer)
palettes <- brewer.pal.info%>%
  mutate(name=row.names(brewer.pal.info))%>%
  filter(category == "qual",colorblind == TRUE)

#vector of 28 colors
col_vector <- unlist(mapply(brewer.pal, palettes$maxcolors, palettes$name))#28 colors

#matching each color with an asv orders present in dataset
col_order <- data.frame(asv_order=unique(zbj_filtered$asv_order), col=col_vector[1:length(unique(zbj_filtered$asv_order))])
head(col_order)

#setting transparency to colors for links 
library(scales)
col_order$col_links <- alpha(col_order$col,0.7)
head(col_order)

#merging col_order with all ASVs 
zbj_col <- unique(merge(zbj_filtered[,c("prey","asv_order"),],col_order,by="asv_order"))
head(zbj_col)

#creating new names for the networks
names(zbj_nets)
new_names <- c("Ayos landscape", "Bokito landscape","Konye landscape")
network_names <- data.frame(net=names(zbj_nets),new_names)
head(network_names)

####plotting a bipartite networks for each landscape
##nets: list of networks to use
##col_asvs: asvs colors to use
##primer: primer that I'm using
for (i in 1:length(zbj_nets)){
  primer <- "zbj"
  net <- zbj_nets[[i]]
  net_name <- network_names[network_names$net ==  names(zbj_nets[i]),]
  orders_plot <- zbj_col[zbj_col$prey %in% rownames(net),]#getting orders from that network
  net_neworder <- net[orders_plot$prey, ] # Reorder rows to match colours by name manually
  png(filename=paste("outputs/plots/network analysis/freq_landscape_network_",primer,"_",net_name$net,".png",sep=""),res= 300, height= 3000, width= 5000)#Open a device, using png()
  plotweb(net_neworder,  ybig=0.15, col.low= orders_plot$col, bor.col.low=orders_plot$col, col.interaction = rep(orders_plot$col_links, each = ncol(net_neworder)), bor.col.interaction=rep(orders_plot$col_links, each = ncol(net_neworder)), labsize=3, adj.high=c(0.5,0), method="normal",low.lablength=0, y.width.low=0.07, y.width.high=0.1)# could also use "cca" for method - avoids overlapping; col.high="orange4" and bor.col.high = "orange4" to change color of  nodes
  text(x = 0.65, y = 1.75, paste(net_name$new_names," (",primer,")",sep=""),cex = 2.5, col = "black")
  legend(x=0,y=0.4,legend = unique(orders_plot$asv_order),x.intersp=0.3, ncol = 4, pch = 15,bty ="n", xpd = TRUE,inset = c(0,0), cex=1.7, title="Orders of arthropod ASVs", text.col = "gray20", title.col = "black",adj = c(0, 0.6),col=unique(orders_plot$col)) 
  dev.off()#clear all plots and saving last plot
}


######Network analyses
###Creating 10000 new networks based on on previous networks for the null model
nullnetworks<-lapply(zbj_nets, nullmodel,N=1000,method="vaznull")#here just using 1000 null networks as an example 
#Method can be changed (see nullmodel guide in bipartite package info)


###Network level metrics
##creating a empty data frame to store all zscores: z-scores are used to know if our metric values are significant compared with the null model
zscores <- data.frame()

library(ggplot2)
##Nestedness - may time some time to compute
#calculating nestedness for each of our networks as weighted nestedness based on overlap and decreasing fill (WNODF). High values indicate more nested networks. 
for (i in 1:length(zbj_nets)){
  net <- zbj_nets[[i]]
  null <- nullnetworks[[i]]
  #measuring nestedness for each network
  nest<-networklevel(net,index="weighted NODF")#
  #measuring nestedness for each of the null models
  nullnest<-unlist(sapply(null,networklevel,index="weighted NODF"))
  #calculates the z-score 
  zscore_nest <- (nest-mean(nullnest))/sd(nullnest)
  #plot to see wheres does our zscore stands compared to the null models
  plot <- ggplot()+
    aes(nullnest)+
    geom_histogram(bins=100,fill="darkolivegreen3",colour="black")+
    theme_classic()+
    geom_vline(xintercept=nest,linetype="solid",linewidth=1,colour="red")+
    ggtitle(names(zbj_nets[i])) 
  print(plot)
  ggsave(paste("outputs/plots/network analysis/null_model_",names(zscore_nest),"_",names(zbj_nets[i]),".png",sep = ""), plot = plot, width =5 , height = 5) #To save the plot
  zscores <- rbind(zscores,data.frame(metric = names(zscore_nest), network=names(zbj_nets[i]),value = nest, zscore=as.vector(zscore_nest)))
}
#The three networks are less nested than expected, with Ayos and Konye being less nested than Bokito 

##Modularity - may time some time to compute
#calculating modularity for each of our networks. High values indicate higher modularity.
for (i in 1:length(zbj_nets)){
  net <- zbj_nets[[i]]
  null <- nullnetworks[[i]]
  #measuring Modularity for each network
  mod<-computeModules(net,method="Beckett",steps=100)#higher number of steps to increase precision - here just using 100 as an example. 1000 steps are advised. 
  #measuring Modularity for each of the null models. 
  nullmod<-sapply(null,computeModules,method="Beckett",steps=100)#higher number of steps to increase precision - here just using 100 as an example. 1000 steps are advised. 
  #Creating object with all likelihoods 
  nulllikeli <- sapply(nullmod,function(x)x@likelihood)
  #calculates the z-score based on the likelihood 
  zscore_mod <- (mod@likelihood-mean(nulllikeli))/sd(nulllikeli)
  #plot to see wheres does our zscore stands compared to the null models
  plot <- ggplot()+
    aes(nulllikeli)+
    geom_histogram(binwidth = 0.001,fill="darkolivegreen3",colour="black")+
    theme_classic()+
    geom_vline(xintercept=mod@likelihood,linetype="solid",linewidth=1,colour="red")+
    ggtitle(names(zbj_nets[i])) 
  print(plot)
  ggsave(paste("outputs/plots/network analysis/null_model_modularity","_",names(zbj_nets[i]),".png",sep = ""), plot = plot, width =5 , height = 5) #To save the plot
  zscores <- rbind(zscores,data.frame(metric = "Modularity", network=names(zbj_nets[i]),value = mod@likelihood, zscore=as.vector(zscore_mod)))
}
#The three networks were more modular than expected, however Ayos was less modular when compared with the other two. 

#In the case you want to plot the modules - not very useful for big networks
mod_Bokito<-computeModules(zbj_nets$Bokito)#example for Bokito landscape
plotModuleWeb(mod_Bokito)

##saving z-scores to csv file
write.csv2(zscores,"./outputs/results/network analysis/network_zscores_zbj.csv",row.names = FALSE)

#plotting s-scores 
metrics_plot <- ggplot(data=zscores, aes(x=network, y=zscore))+ 
  geom_bar(stat='identity')+
  #scale_shape(solid=F)+ #no fill for points
  facet_wrap(~metric,ncol=2,scales="free_y")+
  theme_bw()+ 
  xlab("")+ ylab("Value")+
  theme(axis.title.x=element_text(size=25))+
  theme(axis.title.y=element_text(size=25))+
  theme(axis.text.x = element_text(colour = "black", size= 15, angle = 45,vjust = 1, hjust = 1))+
  theme(axis.text.y=element_text(size=25))+
  theme(legend.text=element_text(size=25), legend.title = element_text(size=25), legend.position="bottom", legend.key = element_blank())+
  theme(strip.text.x = element_text(face="italic",size=25))+
  theme(strip.background = element_rect(fill = 'grey'))
metrics_plot
#We see that Modularity and Nestedness were significant in all networks. Bokito and Konye were more modular while Ayos and Konye were more nestedness. 

ggsave("./outputs/plots/network analysis/zscore_metrics_zbj.png", plot = metrics_plot, width =15 , height = 12) #To save the plot
