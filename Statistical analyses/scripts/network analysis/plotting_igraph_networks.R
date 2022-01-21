##################################################
## Project: metabarcoding and networks
## Script purpose: creating and plotting with igraph 
## Date: 10/09/2021
## Author: Diogo F. Ferreira (ferreiradfa@gmail.com)
## Packages: igraph, qgraph and RColorBrewer
## Notes: you will need the function to generate_final_data.r to run this script
###################################################

###importing metabarcoding data 
#zbj primer: 224 samples, 2 controls and 2889 unique OTUS
zbj <- read.csv2("data/16Farms_OTU_Table_vsearch.csv", dec = ".", row.names = 1,header = TRUE, check.names=FALSE,na.strings=c("NA", "NULL", "", ".")) 
dim(zbj)
colnames(zbj)
head(zbj)

###importing species information for each sample
library(openxlsx)#read excel and sheet
samples_list <- read.xlsx("data/faeces_sample_database.xlsx", sheet = "samples")
colnames(samples_list)
head(samples_list)

###Reading function that returns final data frame with all data organized and cleaned - see organise&clean_metabarcoding.r for more details
source('scripts/organize&clean_metabarcoding.r')

unique(zbj$class)
#remove all OTUs that represent less than 1%, keeping onky Arachnida and Insecta, and remove OTUs not identified until order
zbj_data <- final_metbar(data = zbj,sample_list = samples_list, remove_samples=F,otus_clean=1, keep_class=c("Arachnida","Insecta"),remove_NAorders=T,remove_NAfamily=F,desired_species=NULL)
head(zbj_data)
dim(zbj_data)
#1499 rows and 18 columns
n_distinct(zbj_data$prey)#number of different otus
#816
n_distinct(zbj_data$predator)#number of samples 
#214


####aggregate data by site and species
colnames(zbj_data)#keeping only species, prey, otu_order, landscape, weight, proportion and creating column for otu frequency
zbj_aggr <- zbj_data %>%
  group_by(predator=species_code,prey,otu_order,landscape)%>%
  summarise(weight=sum(weight),proportion=mean(proportion),freq_otus=n())
head(zbj_aggr)
dim(zbj_aggr)
#1214 rows and 7 columns
levels(as.factor(zbj_aggr$predator))#check how and what bats/bird species I have
#11 

####Creating data frame (nodes list) with frequency and abundance of all bat and bird species plus otus 
###edge list
##predators - insects
##how many times each otu occurred per predator species
pred_ins <- zbj_aggr %>%
  group_by(Source=predator,Target=prey,landscape)%>%
  summarise(Weight=sum(freq_otus))
head(pred_ins)

##checking how many samples per species
samples_sp <- zbj_data %>%
  group_by("Source"=species_code,landscape)%>%
  summarise(samples = n_distinct(predator))
head(samples_sp)

##standardizing samples by bat species where I have that otu
pred_ins_std <- merge(pred_ins, samples_sp, by = c("Source","landscape"))
head(pred_ins_std)

pred_ins_std$Weight_std <- pred_ins_std$Weight/pred_ins_std$samples*100
head(pred_ins_std)

edge_final <- pred_ins_std%>%
  select(Source,Target,"weight"=Weight_std,landscape)
head(edge_final)

###node list
##importing capture data
captures <- read.csv("data/bird_and_bat_capturedata.csv",sep = ";")
head(captures)

#captures per landscape
all_abun <- captures %>%
  group_by(id=code,taxon,landscape) %>%
  summarize(count = n())
head(all_abun)

###calculating otus frequency per landscape
#only makes sense if I have the same number of samples per site (not the case, at least for now). It will always be correlated with number of samples per site
otus_freq <- zbj_aggr %>%
  group_by(id=prey,taxon=otu_order,landscape)%>%
  summarise(count=sum(freq_otus))#ideally I should have data for the abundance of insects (check Cyril data in the future)
head(otus_freq)

###join predators abundance and otus frequency
nodes_count <- rbind(all_abun,otus_freq)
head(nodes_count)
dim(nodes_count)

###create a new column in the node list to tell which nodes are the preys and the predators
nodes_count$type <- ifelse(grepl("OTU" , nodes_count$id),"prey","predator")# column in the node list to tell which animal it is
head(nodes_count)
dim(nodes_count)

###selecting only nodes that are present in my samples (just in case we capture species but do not have samples)
nodes_network <- c(unique(zbj_aggr$predator), unique(zbj_aggr$prey))
nodes_final <- nodes_count[nodes_count$id %in% nodes_network,]
head(nodes_final)
dim(nodes_final)

####To generate a list with a igraph network for each site and then for plotting it
library(igraph)
library(qgraph)

igraph_list <- list()#creating empty list

landscape <- unique(zbj_aggr$landscape)# names of all landscape
print(landscape)

for(i in 1:length(landscape)){
  edge_list <- edge_final[which(edge_final$landscape==landscape[i]),]
  nodes <- unique(data.frame(id = c(pull(edge_list,colnames(edge_list[1])),pull(edge_list,colnames(edge_list[2])))))###nodes list contains all unique otus and predators ID
  nodes_info <- nodes_final[which(nodes_final$landscape==landscape[i]),]
  nodes_list <- merge(nodes, nodes_info, by=c("id"))#merging nodes with nodes_count;
  nodes_list <- nodes_list[order(nodes_list$type),]
  network <- graph_from_data_frame(d=edge_list, vertices=nodes_list, directed=F)#creating an igraph object based on the node and edge lists; In igraph a bipartite network is one that has a type vertex attribute; The 'B' letter after "UNW" tells you that this is a bipartite graph
  igraph_list[[i]] <- network
}

names(igraph_list) <- landscape #matching data.frames names to landscape
names(igraph_list)

##Exploring igraph - doing an example for Bokito landscape
nets_igraph <- igraph_list$Bokito
class(nets_igraph)
nets_igraph

#Count number of edges
gsize(nets_igraph)
#427

#Count number of vertices
gorder(nets_igraph)
#337

#View attributes of first five vertices in a dataframe
V(nets_igraph)[[1:5]] 

#View attributes of first five edges in a dataframe
E(nets_igraph)[[1:5]] 

##plotting 
#Set node shapes by type
V(nets_igraph)$shape <- V(nets_igraph)$type
V(nets_igraph)$shape <- gsub("predator","circle",V(nets_igraph)$shape)
V(nets_igraph)$shape <- gsub("prey","square",V(nets_igraph)$shape)

#if we had another variable we could change colors of edges - for example by season
#E(nets_igraph)$color <- E(nets_igraph)$season
#E(nets_igraph)$color <- gsub("wet","#D05B9B",E(nets_igraph)$color)
#E(nets_igraph)$color <- gsub("dry","#5770AE",E(nets_igraph)$color)

#Set nodes colours by taxon
library("RColorBrewer")
taxa <- unique(V(nets_igraph)$taxon)#insect orders plus predator 
col <- brewer.pal(n = length(taxa), name = "RdBu")#one different color for each order
V(nets_igraph)$color <- V(nets_igraph)$taxon

#renaming taxon
V(nets_igraph)$taxon <- gsub("bat","Bat", V(nets_igraph)$taxon)


#replacing order by the specific colour
for(i in 1:length(taxa)){
  V(nets_igraph)$color <- gsub(taxa[i],col[i],V(nets_igraph)$color)
}

###creating layouts for plot
l <- layout_nicely(nets_igraph)#a common layout for nodes
#or
e <- get.edgelist(nets_igraph,names=FALSE)
lnodes <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(nets_igraph))#avoids overlapping nodes

#avoids overlapping edges by curving them
curvesG <- curve_multiple(nets_igraph)

###Plotting igraph object
png(filename="outputs/plots/igraph_network_bokito.png",res= 300, height= 3000, width= 3500)

plot(nets_igraph, 
     vertex.color = V(nets_igraph)$color, 
     vertex.frame.color= V(nets_igraph)$color, 
     vertex.shape = V(nets_igraph)$shape,
     vertex.size=10,
     vertex.label=V(nets_igraph)$taxon,
     vertex.label.color="black",
     vertex.label.cex=.5,
     edge.color = E(nets_igraph)$color,
     edge.width = log(E(nets_igraph)$weight),
     edge.curved=curvesG,
     layout=lnodes)#I can also use "l" 

##I need to redo this - Create a legend for the plot 
legend(x = 1,
       y = 0.5, 
       legend = unique(V(nets_igraph)$taxon),
       fill = unique(V(nets_igraph)$color), 
       title="Taxon",
       text.col = "gray20", 
       title.col = "black", 
       box.lwd = 0, 
       cex = 1,
       bty ="n")
title(main = "Arthropd-Predator network - Bokito", cex.main=2)

dev.off()#clear all plots and saving last plot
