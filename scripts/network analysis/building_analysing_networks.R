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

zbj_data <- final_metbar(data = zbj,sample_list = samples_list, remove_samples=F,otus_clean=T,keep_class=NULL,remove_NAorders=T,remove_NAfamily=F,desired_species=NULL)
head(zbj_data)
dim(zbj_data)
#403 rows and 18 columns
n_distinct(zbj_data$prey)#number of different otus
#196
n_distinct(zbj_data$predator)#number of samples - originally 224 samples 
#166


####aggregate data by site and species
colnames(zbj_data)#keeping only species, prey, otu_order, landscape, weight, proportion and creating column for otu frequency
zbj_aggr <- zbj_data %>%
    group_by(predator=species_code,prey,otu_order,landscape)%>%
    summarise(weight=sum(weight),proportion=mean(proportion),freq_otus=n())
head(zbj_aggr)
dim(zbj_aggr)
#470 rows and 7 columns
levels(as.factor(zbj_aggr$predator))#check how and what bats/bird species I have
#11 

##checking how many samples per species and landscape
zbj_species <- zbj_data %>%
    group_by(landscape,animal,species_code)%>%
    summarise(samples = n_distinct(predator))%>%
    mutate(total_samples = sum(samples),n_species = n_distinct(species_code))%>%#creates columns with total samples per landscaoe and total number of species
    rename("predator"=species_code)
head(zbj_species)

write.csv(zbj_species,"outputs/data/network analysis/n_samples_network_zbj.csv",row.names = FALSE)

##I'll remove species that have less than 2 samples (ideally the cut-off should be higher - at least 10 samples) per landscape
species_remove_zbj <- unique(zbj_species[zbj_species$samples < 10, ])
levels(as.factor(species_remove_zbj$predator))
#I need to remove 3 species

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
####more specific analysis
###Creating 10000 new networks based on our networks for the null model
#should I create null models for every test or can I always use the same?
nullnetworks<-lapply(zbj_nets, nullmodel,N=10000,method="r2dtable")
#r2dtable: Generates null matrices based on Patefield's method.
#vaznull: Generates null matrices based on the Vazquez's method (connectance is kept constant).
#Method can be changed (see nullmodel guide in bipartite), but I personally find Vazquez's method tends to provide the most meaningful comparison.

#I can also use CI at 95%:
Vazquez.null.wNODF<-unlist(sapply(zbj_nets, networklevel, index="weighted NODF"))
mean(Vazquez.null.wNODF)
quantile(unlist(Vazquez.null.wNODF), c(0.025, 0.975))


###Network/Group level metrics
##creating a empty data frame to store all zscores: z-scores are used to know if our metric values are significant compared with the null model
zscores <- data.frame()

##Generality 
#Weighted mean effective number of low level species per high level species (generality), weighted by their marginal totals (row sums). See Tylianakis et al. (2007) and Bersier et al. (2002). This is identical to exp(“partner diversity”, i.e., simply the Jost (2006)-recommended version of diversity
#lower values mean that networks are highly specialized
library(ggplot2)
for (i in 1:length(zbj_nets)){
  net <- zbj_nets[[i]]
  null <- nullnetworks[[i]]
  #measuring generality for each network
  gener<-grouplevel(net,index="generality",weighted=TRUE,level="higher")
  #measuring generality for each of the null models
  nullgener<-unlist(sapply(null,grouplevel,index="generality",weighted=TRUE,level="higher"))
  #calculates the z-score
  zscore_gener <- (gener-mean(nullgener))/sd(nullgener) ##z-score is negative - what does it means? Should I create binary networks because of that?
  #plot to see wheres does our zscore stands compared to the null models
  plot <- ggplot()+
    aes(nullgener)+
    geom_histogram(bins=100,colour="coral2",fill="darkolivegreen3")+
    theme_classic()+
    geom_vline(xintercept=gener,linetype="dashed",size=0.1,colour="red")+
    ggtitle(names(zbj_nets[i])) 
  print(plot)
  ggsave(paste("outputs/plots/network analysis/generality/null_model_",names(zscore_gener),"_",names(zbj_nets[i]),".png",sep = ""), plot = plot, width =5 , height = 5) #To save the plot
  zscores <- rbind(zscores,data.frame(metric = names(zscore_gener), network=names(zbj_nets[i]),value = gener, zscore=as.vector(zscore_gener)))
}


##Vulnerability
#Weighted mean effective number of HL species per LL species (vulnerability), weighted by their marginal totals (row sums). See Tylianakis et al. (2007) and Bersier et al. (2002). This is identical to exp(“partner diversity”, i.e., simply the Jost (2006)-recommended version of diversity. 
#Lower values mean that networks are highly specialized
for (i in 1:length(zbj_nets)){
  net <- zbj_nets[[i]]
  null <- nullnetworks[[i]]
  #measuring vulnerability for each network
  vuln<-grouplevel(net,index="vulnerability",weighted=TRUE,level="lower")#
  #measuring vulnerability for each of the null models
  nullvuln<-unlist(sapply(null,grouplevel,index="vulnerability",weighted=TRUE,level="lower"))
  #calculates the z-score 
  zscore_vuln <- (vuln-mean(nullvuln))/sd(nullvuln)
  #plot to see wheres does our zscore stands compared to the null models
  plot <- ggplot()+
    aes(nullvuln)+
    geom_histogram(bins=100,colour="coral2",fill="darkolivegreen3")+
    theme_classic()+
    geom_vline(xintercept=vuln,linetype="dashed",size=0.1,colour="red")+
    ggtitle(names(zbj_nets[i])) 
  print(plot)
  ggsave(paste("outputs/plots/network analysis/vulnerability/null_model_",names(zscore_vuln),"_",names(zbj_nets[i]),".png",sep = ""), plot = plot, width =5 , height = 5) #To save the plot
  zscores <- rbind(zscores,data.frame(metric = names(zscore_vuln), network=names(zbj_nets[i]),value = vuln, zscore=as.vector(zscore_vuln)))
}


##Interaction evenness
#Shannon’s evenness for the web entries. Note that the two options are rather different. By definition, IE = H/Hmax; H = -sum(p.i.mat*log(p.i.mat)), where p.i.mat = matrix/sum(entries in matrix). This means, when calculating H, do we treat all possible links as species, and the interactions (cell values) as measure of their abundance? By definition, Hmax = ln(N). The key question is: What is the right value for N? Since we treat the matrix cells as species, it is (clearly?) the number of matrix cells, i.e. number of higher trophic level species x number of lower trophic level species. We think this logic justifies our default "prod". However, others argue in favour of N=number of links. Please see note for our discussion on this point.
#Lower values mean that networks are highly specialized
for (i in 1:length(zbj_nets)){
  net <- zbj_nets[[i]]
  null <- nullnetworks[[i]]
  #measuring interaction evenness for each network
  even<-networklevel(net,index="interaction evenness",weighted=TRUE)
  #measuring interaction evenness for each of the null models
  nulleven<-unlist(sapply(null,networklevel,index="interaction evenness",weighted=TRUE))
  #calculates the z-score
  zscore_even <- (even-mean(nulleven))/sd(nulleven)
  #plot to see wheres does our zscore stands compared to the null models
  plot <- ggplot()+
    aes(nulleven)+
    geom_histogram(bins = 100,colour="coral2",fill="darkolivegreen3")+
    theme_classic()+
    geom_vline(xintercept=even,linetype="dashed",size=0.1,colour="red")+
    ggtitle(names(zbj_nets[i])) 
  print(plot)
  ggsave(paste("outputs/plots/network analysis/interaction evenness/null_model_",names(zscore_even),"_",names(zbj_nets[i]),".png",sep = ""), plot = plot, width =5 , height = 5) #To save the plot
  zscores <- rbind(zscores,data.frame(metric = names(zscore_even), network=names(zbj_nets[i]),value = even,zscore=as.vector(zscore_even)))
}

##Specialization H2'
#H2’ is an index describing the level of “complementarity specialisation” (or should one say: selectiveness?) of an entire bipartite network (Blüthgen et al. 2006). It describes to which extent observed interactions deviate from those that would be expected given the species marginal totals.
#Higher values mean that networks are highly specialized, i.e., more selective a species the larger is H2’ for the web
for (i in 1:length(zbj_nets)){
  net <- zbj_nets[[i]]
  null <- nullnetworks[[i]]
  #measuring specialization for each network
  h2<-H2fun(net,H2_integer=TRUE)
  #measuring specialization for each of the null models
  nullh2<-sapply(null,H2fun,H2_integer=TRUE)
  #calculates the z-score
  zscore_h2 <- (h2[1]-mean(nullh2[1,]))/sd(nullh2[1,])
  #plot to see wheres does our zscore stands compared to the null models
  plot <- ggplot()+
    aes(nullh2[1,])+
    geom_histogram(bins = 100,colour="coral2",fill="darkolivegreen3")+
    theme_classic()+
    geom_vline(xintercept=h2[1],linetype="dashed",size=1,colour="red")+
    ggtitle(names(zbj_nets[i])) 
  print(plot)
  ggsave(paste("outputs/plots/network analysis/h2/null_model_",names(zscore_h2),"_",names(zbj_nets[i]),".png",sep = ""), plot = plot, width =5 , height = 5) #To save the plot
  zscores <- rbind(zscores,data.frame(metric = names(zscore_h2), network=names(zbj_nets[i]),value = h2[1],zscore=as.vector(zscore_h2)))
}

##Modularity
#calculating modularity for each of our networks - using a higher number of stepps to inrease precision
for (i in 1:length(zbj_nets)){
  net <- zbj_nets[[i]]
  null <- nullnetworks[[i]]
  #measuring Modularity for each network
  mod<-computeModules(net,method="Beckett",steps=1000)
  #measuring Modularity for each of the null models
  nullmod<-sapply(null,computeModules,method="Beckett",steps=1000)
  #Creating object with all likelihood 
  nulllikeli <- sapply(nullmod,function(x)x@likelihood)
  #calculates the z-score based on the likelihood 
  zscore_mod <- (mod@likelihood-mean(nulllikeli))/sd(nulllikeli)##z-score is negative - what does it means? Should I create binary networks because of that?
  #plot to see wheres does our zscore stands compared to the null models
  plot <- ggplot()+
    aes(nulllikeli)+
    geom_histogram(binwidth = 0.001,colour="coral2",fill="darkolivegreen3")+
    theme_classic()+
    geom_vline(xintercept=mod@likelihood,linetype="dashed",size=1,colour="red")+
    ggtitle(names(zbj_nets[i])) 
  print(plot)
  ggsave(paste("outputs/plots/network analysis/modularity/null_model_",names(zscore_mod),"_",names(zbj_nets[i]),".png",sep = ""), plot = plot, width =5 , height = 5) #To save the plot
  zscores <- rbind(zscores,data.frame(metric = "Modularity", network=names(zbj_nets[i]),value = mod@likelihood, zscore=as.vector(zscore_mod)))
}


#in the case I want to plot the modules 
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

###Species level metrics
##Specialization D'
#The d’ index is derived from Kulback-Leibler distance (as is Shannon’s diversity index), and calculates how strongly a species deviates from a random sampling of interacting partners available. It ranges from 0 (no specialization) to 1 (perfect specialist). In the case of a pollination web, a pollinator may be occurring only on one plant species, but if this species is the most dominant one, there is limited evidence for specialization. Hence this pollinator would receive a low value. In contrast,a pollinator that occurs only on the two rarest plants would have a very high value of d’.
#The way this function is implemented, it calculates expected values for each cell based on the product of observed marginal sums (i.e. column and row sums) times sum(web). Then it rounds off to integers and allocates the remaining interactions in two steps: First, all columns and rows with marginal sums of 0 obtain one interaction into the cell with the highest expected value. Secondly, all remaining interactions are distributed according to difference between present and expected value: those cells with highest discrepancy receive an interaction until the sum of all entries in the new web equals those in the original web. Now the d-values for this web are calculated and used as dmin.
#dfun returns the d’ values for the lower trophic level. Use fun(t(web)) to get the d’-values for the higher trophic level (as does specieslevel). If you want to provide external abundances, you must provide those of the OTHER trophic level! (This help file is written as if you were interested in the lower trophic level.)
#Higher values mean that species are highly specialized

d_species <- data.frame()#to store d' values
for (i in 1:length(zbj_nets)){
  net <- t(zbj_nets[[i]])#i need to transpose to get d' for my higher level (birds)
  #measuring specialization for each network
  d<-data.frame(dfun(net, abuns=NULL))
  d$species <- row.names(d)
  #save file
  d_species <- rbind(d_species,data.frame(d,network=names(zbj_nets[i])))
}

#changing order of columns
colnames(d_species)
d_species <- d_species[c(6,5,1:4)]#reordering data frame

##saving d value to csv file
write.csv2(d_species,"outputs/results/network analysis/network_dvalues_zbj.csv",row.names = FALSE)


##plotting d values
library(ggthemes)
#The errorbars overlapped, so use position_dodge to move them horizontally
plot_d <- ggplot(data=d_species, aes(x=network, y=dprime))+ #again, excluding intercept because estimates so much larger
  geom_bar(stat='identity')+
  #scale_shape(solid=F)+ #no fill for points
  facet_wrap(~species,ncol=2)+#,scales="free_y")+
  theme_bw()+ 
  xlab("")+ ylab("d'bird")+
  theme(axis.title.x=element_text(size=25))+
  theme(axis.title.y=element_text(size=25))+
  theme(axis.text.x = element_text(colour = "black", size= 15, angle = 45,vjust = 1, hjust = 1))+
  theme(axis.text.y=element_text(size=25))+
  theme(legend.text=element_text(size=25), legend.title = element_text(size=25), legend.position="bottom", legend.key = element_blank())+
  theme(strip.text.x = element_text(face="italic",size=25))+
  theme(strip.background = element_rect(fill = 'grey'))
plot_d

ggsave("outputs/plots/network analysis/network_dvalues_zbj.png", plot = plot_d, width =15 , height = 15) #To save the plot

