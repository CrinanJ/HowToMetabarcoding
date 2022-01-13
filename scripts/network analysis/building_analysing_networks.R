##################################################
## Project: metabarcoding and networks
## Script purpose: creating, plotting and analyzing bipartite networks 
## Date: 10/09/2021
## Author: Diogo F. Ferreira (ferreiradfa@gmail.com)
## Packages: bipartite, ggplot2, ggthemes, GISTools and RColorBrewer 
## Notes: you will need the function to generate_final_data.r to run this script
###################################################

###importing metabarcoding data 
#zbj primer: 224 samples and 3105 unique OTUS
zbj <- read.csv2("data/16Farms_OTU_Table_withConfidence.csv", dec = ".", row.names = 1,header = TRUE, check.names=FALSE,na.strings=c("NA", "NULL", "", ".")) 
colnames(zbj)
head(zbj)

###importing species information for each sample
library(openxlsx)#read excel and sheet
samples_list <- read.xlsx("data/faeces_sample_database.xlsx", sheet = "samples")
colnames(samples_list)
head(samples_list)

###Reading function that returns final data frame with all data organized and cleaned - see organise&clean_metabarcoding.r for more details
source('scripts/organize&clean_metabarcoding.r')

zbj_data <- final_metbar(data = zbj,sample_list = samples_list, remove_samples=F,otus_clean=0, keep_class=NULL,remove_NAorders=F,remove_NAfamily=F,desired_species=NULL)
head(zbj_data)
dim(zbj_data)
#6426 rows and 18 columns
n_distinct(zbj_data$prey)#number of different otus
#3099
n_distinct(zbj_data$predator)#number of samples - originally 224 samples 
#220


####aggregate data by site and species
colnames(zbj_data)#keeping only species, prey, otu_order, landscape, weight, proportion and creating column for otu frequency
zbj_aggr <- zbj_data %>%
    group_by(predator=species_code,prey,otu_order,landscape)%>%
    summarise(weight=sum(weight),proportion=mean(proportion),freq_otus=n())
head(zbj_aggr)
dim(zbj_aggr)
#319 rows and 7 columns
levels(as.factor(zbj_aggr$predator))#check how and what bats/bird species I have
#11 

##checking how many samples per species and landscape
zbj_species <- zbj_data %>%
    group_by(landscape,animal,species_code)%>%
    summarise(samples = n_distinct(predator))%>%
    mutate(total_samples = sum(samples),n_species = n_distinct(species_code))%>%#creates columns with total samples per landscaoe and total number of species
    rename("predator"=species_code)
head(zbj_species)

write.csv(zbj_species,"outputs/results/network analysis/n_samples_network_zbj.csv",row.names = FALSE)

##I'll remove species that have less than 5 samples (ideally higher the cut-off the better) per landscape 
species_remove_zbj <- unique(zbj_species[zbj_species$samples < 5, ])
levels(as.factor(species_remove_zbj$predator))
#I need to remove 5 species

zbj_filtered <- data.frame()
for (i in 1:length(unique(zbj_aggr$landscape))){
  land <- zbj_aggr[zbj_aggr$landscape==unique(zbj_aggr$landscape)[i],]#select rows from a landscape
  remove <- species_remove_zbj[species_remove_zbj$landscape == unique(zbj_aggr$landscape)[i],]$predator
  final <- land[!land$predator %in% remove,] #removes species with only one sample
  zbj_filtered <- rbind(zbj_filtered, final)
}
head(zbj_filtered)
dim(zbj_filtered)
#339 rows and 7 columns
levels(as.factor(zbj_filtered$predator))
#9 species

########################Bipartite networkt - on bipartite package#################
####generating an network based on a adjacency matrix and creating an each network based on webID /farm for now)
#if I want to generate networks for each site and season, I need to change webID and add "paste(zbj_filtered$farm,zbj_filtered$season,sep="-")" to webID
library(bipartite)
zbj_nets <- frame2webs(data.frame(lower=zbj_filtered$prey,higher=zbj_filtered$predator,webID=zbj_filtered$landscape,freq=zbj_filtered$freq_otus),type.out="list")
lapply(zbj_nets, dim)#checking dimensions of all data frames

unlist(lapply(zbj_nets, function(x) sum(x>=1)))#number of links per network
#Ayos Bokito  Konye 
#182    121    154 

unlist(lapply(zbj_nets, sum))#total otus frequency per network
#Ayos Bokito  Konye 
#191    202    179 

unlist(lapply(zbj_nets, nrow))#number of otus per network
#Ayos Bokito  Konye 
#162    106    122 

unlist(lapply(zbj_nets, ncol))#number of species per network
#Ayos Bokito  Konye 
#7      4      8 

#####This tiny chunk shows how many MOTU are shared between networks
samples <- lapply(zbj_nets, rownames)
Reduce(intersect, samples)
length(Reduce(intersect, samples))
#18 otus shared
(length(Reduce(intersect, samples))/length(unique(zbj_filtered$prey)))*100
#6.020067% of otus are shared between all networks

####Plotting bipartite network
##creating pallete (28 colours and colorblind) to give to each arthrpod order
library(RColorBrewer)
palettes <- brewer.pal.info%>%
  mutate(name=row.names(brewer.pal.info))%>%
  filter(category == "qual",colorblind == TRUE)

#vector of 28 colours
col_vector <- unlist(mapply(brewer.pal, palettes$maxcolors, palettes$name))#28 colours

#matching each colour with an otu orders present in dataset
col_order <- data.frame(otu_order=unique(zbj_filtered$otu_order), col=col_vector[1:length(unique(zbj_filtered$otu_order))])
head(col_order)

#setting transparency to colours for links
library(GISTools)
col_order$col_links <- add.alpha(col_order$col,0.7)
head(col_order)

#mergind col_order with all OTUs 
zbj_col <- unique(merge(zbj_filtered[,c("prey","otu_order"),],col_order,by="otu_order"))
head(zbj_col)

#creating new names for the networks
names(zbj_nets)
new_names <- c("Ayos landscape", "Bokito landscape","Konye landscape")
network_names <- data.frame(net=names(zbj_nets),new_names)
head(network_names)

####plotting a bipartite networks for each landscape
##nets: list of networks to use
##col_otus: otus colours to use
##primer: primer that I'm using
for (i in 1:length(zbj_nets)){
  primer <- "zbj"
  net <- zbj_nets[[i]]
  net_name <- network_names[network_names$net ==  names(zbj_nets[i]),]
  orders_plot <- zbj_col[zbj_col$prey %in% rownames(net),]#getting orders from that network
  png(filename=paste("outputs/plots/network analysis/freq_landscape_network_",primer,"_",net_name$net,".png",sep=""),res= 300, height= 3000, width= 5000)#Open a device, using png()
  plotweb(net,  ybig=0.15, col.low= orders_plot$col, bor.col.low=orders_plot$col, col.interaction = rep(orders_plot$col_links, each = ncol(net)), bor.col.interaction=rep(orders_plot$col_links, each = ncol(net)), labsize=3, adj.high=c(0.5,0), method="normal",low.lablength=0, y.width.low=0.07, y.width.high=0.1)# could also use "cca" for method - avoids overlapping; col.high="orange4" and bor.col.high = "orange4" to change color of  nodes
  text(x = 0.65, y = 1.75, paste(net_name$new_names," (",primer,")",sep=""),cex = 2.5, col = "black")
  legend(x=0,y=0.4,legend = unique(orders_plot$otu_order),x.intersp=0.3, ncol = 4, pch = 15,bty ="n", xpd = TRUE,inset = c(0,0), cex=1.7, title="Orders of arthropod OTUs", text.col = "gray20", title.col = "black",adj = c(0, 0.6),col=unique(orders_plot$col)) 
  dev.off()#clear all plots and saving last plot
}


######Network analysis
###Creating 10000 new networks based on on previous networks for the null model
nullnetworks<-lapply(zbj_nets, nullmodel,N=1000,method="vaznull")
#Method can be changed (see nullmodel guide in bipartite package info)

###Network level metrics
##creating a empty data frame to store all zscores: z-scores are used to know if our metric values are significant compared with the null model
zscores <- data.frame()

library(ggplot2)
##Nestedness
#calculating nestedness for each of our networks as weighted nestedness based on overlap and decreasing fill (WNODF). High values indicate nestedness.
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
    geom_vline(xintercept=nest,linetype="solid",size=1,colour="red")+
    ggtitle(names(zbj_nets[i])) 
  print(plot)
  ggsave(paste("outputs/plots/network analysis/null_model_",names(zscore_nest),"_",names(zbj_nets[i]),".png",sep = ""), plot = plot, width =5 , height = 5) #To save the plot
  zscores <- rbind(zscores,data.frame(metric = names(zscore_nest), network=names(zbj_nets[i]),value = nest, zscore=as.vector(zscore_nest)))
}


##Modularity
#calculating modularity for each of our networks. High values indicate modularity.
for (i in 1:length(zbj_nets)){
  net <- zbj_nets[[i]]
  null <- nullnetworks[[i]]
  #measuring Modularity for each network
  mod<-computeModules(net,method="Beckett",steps=1000)#higher number of steps to increase precision
  #measuring Modularity for each of the null models
  nullmod<-sapply(null,computeModules,method="Beckett",steps=1000)#higher number of steps to increase precision
  #Creating object with all likelihood 
  nulllikeli <- sapply(nullmod,function(x)x@likelihood)
  #calculates the z-score based on the likelihood 
  zscore_mod <- (mod@likelihood-mean(nulllikeli))/sd(nulllikeli)##z-score is negative - what does it means? Should I create binary networks because of that?
  #plot to see wheres does our zscore stands compared to the null models
  plot <- ggplot()+
    aes(nulllikeli)+
    geom_histogram(binwidth = 0.001,fill="darkolivegreen3",colour="black")+
    theme_classic()+
    geom_vline(xintercept=mod@likelihood,linetype="solid",size=1,colour="red")+
    ggtitle(names(zbj_nets[i])) 
  print(plot)
  ggsave(paste("outputs/plots/network analysis/null_model_modularity","_",names(zbj_nets[i]),".png",sep = ""), plot = plot, width =5 , height = 5) #To save the plot
  zscores <- rbind(zscores,data.frame(metric = "Modularity", network=names(zbj_nets[i]),value = mod@likelihood, zscore=as.vector(zscore_mod)))
}

#In the case you want to plot the modules - not very useful for big networks
mod_Bokito<-computeModules(zbj_nets$Bokito)
plotModuleWeb(mod_Bokito)

##saving z-scores to csv file
write.csv2(zscores,"outputs/results/network analysis/network_zscores_zbj.csv",row.names = FALSE)

metrics_plot <- ggplot(data=zscores, aes(x=network, y=value))+ #again, excluding intercept because estimates so much larger
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

ggsave("outputs/plots/network analysis/network_metrics_zbj.png", plot = metrics_plot, width =15 , height = 15) #To save the plot
