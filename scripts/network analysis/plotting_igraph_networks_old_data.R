##################################################
## Project: bats and networks
## Script purpose: plotting igraph network - igraph package
## Date: 05/01/2021
## Author: Diogo F. Ferreira (ferreiradfa@gmail.com)
## Notes:
###################################################

###importing metabarcoding data 
#zbj primer - most otus of targeted classes (14301 out of 14927) - I have to remove 49 samples due to controls - 269 samples - For Insecta and Arachnida after cleaning and removing: 309 otus to order of 895, and 232 to family.
zbj <- read.csv("data/metabarcoding_old_data/MiSeq1&MiSeq2_16 farms_13April2020/zbj_consensus.tsv", sep = "\t",row.names = 1,header = TRUE, na.strings=c("NA", "NULL", "", ".")) 

#####using zbj - more data overall
###Reading my function that returns final data frame with all data organised and cleaned - see generate_final_data.r for more details
source('scripts/generate_final_data.r')

met_data <- function(met){
  final_metbar(data = met,remove_samples=T,otus_clean=T,keep_class=c("Arachnida","Insecta"),remove_NAorders=T,remove_NAfamily=F,desired_species=NULL)
}


zbj_data <- met_data(zbj)
head(zbj_data)
dim(zbj_data)
n_distinct(zbj_data$prey)#number of different otus
#309
n_distinct(zbj_data$predator)#number of samples
#219

####aggregate data by site and species - keeping only bats and insectivores
met_aggr <- function(met){
  met %>%
    filter(animal == "bat")%>%
    group_by(predator=species,prey,otu_order,pest,landscape,farm)%>%
    summarise(weight=sum(weight),proportion=mean(proportion),freq_otus=n())
}

zbj_aggr <- met_aggr(zbj_data)
head(zbj_aggr)
dim(zbj_aggr)
levels(as.factor(zbj_aggr$predator))
#21 species

####Creating data frame (nodes list) with frequency and abundance of all bat and bird species plus otus 
###importing bird capture data - I am nothing use birds for now but it will be useful for latter
birds <- read.xlsx("data/BirdBandingDataAll_19Oct20.xlsx", sheet = "capdata")
head(birds)

birds$season <- ifelse(birds$month > 6, "wet", "dry")#add column for season 
birds$taxon <- "bird" #new columns to identify type of animal

#bird abundance per site and season
birds_abun <- birds %>%
  group_by(id=species,taxon,landscape,farm) %>%
  summarize(count = n())
head(birds_abun)

###importing bat capture data
bats <- read.xlsx("data/BatCaptureData_Cameroon__20Dec07.xlsx", sheet = "Captures")
head(bats)

#matching column names with birds capture data
names(bats)<-tolower(names(bats))#all column names to lowercase
bats$season<-tolower(bats$season)#season needs to be in lower case
names(bats)[names(bats) == "site"] <- "farm"
names(bats)[names(bats) == "location"] <- "landscape"

bats$taxon <- "bat"#new columns to identify type of animal

#bat abundance per site and season
bats_abun <- bats %>%
  group_by(id=species,taxon,landscape,farm) %>%
  summarize(count = n())
head(bats_abun)

###calculating otus frequency per site and season
#only makes sense if I have the same number of samples per site (not the case, at least for now). It will always be correlated with number of samples per site
otus_freq <- zbj_aggr %>%
  group_by(id=prey,taxon=otu_order,landscape,farm)%>%
  summarise(count=sum(freq_otus))#ideally I should have data for the abundance of insects (check Cyril data in the future)
head(otus_freq)

###join predators abundance and otus frequency
nodes_count <- rbind(birds_abun,bats_abun,otus_freq)
head(nodes_count)
dim(nodes_count)

###create a new column in the node list to tell which nodes are the preys and the predators
nodes_count$type <- ifelse(grepl("OTU" , nodes_count$id),"prey","predator")# column in the node list to tell which animal it is
head(nodes_count)
dim(nodes_count)

###selecting only nodes that are present in my samples - doesn't work super well because I have nodes that only occur in a farm
nodes_network <- c(unique(zbj_aggr$predator), unique(zbj_aggr$prey))
nodes_final <- nodes_count[nodes_count$id %in% nodes_network,]
head(nodes_final)
dim(nodes_final)

write.csv(nodes_final,"outputs/data/gephi_nodelist_zbj.csv",row.names = FALSE)

####creating the edge list
edge_final <- zbj_aggr %>%
  ungroup()%>%
  select(Source=predator,Target=prey,"weight"=freq_otus,landscape,farm)
head(edge_final)

write.csv(edge_final,"outputs/data/gephi_edgelist_zbj.csv",row.names = FALSE)

####To generate a list with a igraph network for each site and then for plotting it
library(igraph)
library(qgraph)

igraph_list <- list()#creating empty list

farms <- unique(zbj_aggr$farm)# names of all farms
print(farms)

for(i in 1:length(farms)){
  edge_list <- edge_final[which(edge_final$farm==farms[i]),]
  nodes <- unique(data.frame(id = c(pull(edge_list,colnames(edge_list[1])),pull(edge_list,colnames(edge_list[2])))))###nodes list contains all unique otus and predators ID
  nodes_info <- nodes_final[which(nodes_final$farm==farms[i]),]
  nodes_list <- merge(nodes, nodes_info, by=c("id"))#merging nodes with nodes_count;
  nodes_list <- nodes_list[order(nodes_list$type),]
  network <- graph_from_data_frame(d=edge_list, vertices=nodes_list, directed=F)#creating an igraph object based on the node and edge lists; In igraph a bipartite network is one that has a type vertex attribute; The 'B' letter after "UNW" tells you that this is a bipartite graph
  igraph_list[[i]] <- network
}

names(igraph_list) <- farms#matching data.frames names to farms
names(igraph_list)

##Exploring igraph - I'll have then to create a loop to apply everything to igraph_lisy
nets_igraph <- igraph_list$NGUI002
class(nets_igraph)
nets_igraph

# Count number of edges
gsize(nets_igraph)
# Count number of vertices
gorder(nets_igraph)
# View attributes of first five vertices in a dataframe
V(nets_igraph)[[1:5]] 
E(nets_igraph)[[1:5]] 

##plotting 
# Set node shapes by type
V(nets_igraph)$shape <- V(nets_igraph)$type
V(nets_igraph)$shape <- gsub("predator","circle",V(nets_igraph)$shape)
V(nets_igraph)$shape <- gsub("prey","square",V(nets_igraph)$shape)

#Set link colors by pest or not pest
#E(nets_igraph)$color <- E(nets_igraph)$pest
#E(nets_igraph)$color <- gsub("no","#D05B9B",E(nets_igraph)$color)
#E(nets_igraph)$color <- gsub("yes","#5770AE",E(nets_igraph)$color)

#Set nodes colours by taxon
library("RColorBrewer")
taxa <- unique(V(nets_igraph)$taxon)#insect orders plus predator 
col <- brewer.pal(n = length(taxa), name = "RdBu")#one different colour for each order
V(nets_igraph)$color <- V(nets_igraph)$taxon

#V(nets_igraph)$taxon <- gsub("bat","Bat", V(nets_igraph)$taxon)
#V(nets_igraph)$taxon <- gsub("bird","Bird",V(nets_igraph)$taxon)

#replacing order by the specific colour
for(i in 1:length(taxa)){
  V(nets_igraph)$color <- gsub(taxa[i],col[i],V(nets_igraph)$color)
}

###creating layouts for plot
#position of nodes for the plot - similar to plotweb
animal <- table(V(nets_igraph)$type)
animal
position <- matrix(c(c(1:animal[1],1:animal[2]),rep(1:2,c(animal[1],animal[2]))),ncol=2,nrow=gorder(nets_igraph))#could also use  sum(animal)

#or I can use the layout layout=layout_as_bipartite() - although doesn't work always
bipartite <- layout.bipartite(nets_igraph)# doesn't work - don't know why

#other layouts
l <- layout_nicely(nets_igraph)#if I want another shape for the plot

e <- get.edgelist(nets_igraph,names=FALSE)
lnodes <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(nets_igraph))#avoids overlapping nodes

#avoids overlapping edges by curving them
curvesG <- curve_multiple(nets_igraph)

###Plotting igraph object
png(filename="outputs/plots/network_ngui2.png",res= 300, height= 3000, width= 3500)

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
     layout=lnodes)#I can also use "l" or "position"

##I need to redo this - Create a legend for the plot 
legend(x = 1.05,
       y = 1, 
       legend = c("Prey","Predator"), 
       pch = c(15,16),  
       title="Taxa", 
       text.col = "gray20", 
       title.col = "black", 
       box.lwd = 0,
       cex = 1.5,
       bty ="n")
legend(x = 1,
       y = 0.65, 
       legend = unique(V(nets_igraph)$taxon),
       fill = unique(V(nets_igraph)$color), 
       title="Taxon",
       text.col = "gray20", 
       title.col = "black", 
       box.lwd = 0, 
       cex = 1.5,
       bty ="n")
title(main = "Arthropd-Predator network - NGUI2", cex.main=2)

dev.off()#clear all plots and savin last plot
