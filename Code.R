#################### Part I #################### 

# Read csv files
df1 <- read.csv2("2009.07.01.csv", sep = ",", header= TRUE)
df2 <- read.csv2("2009.07.02.csv", sep = ",", header= TRUE)
df3 <- read.csv2("2009.07.03.csv", sep = ",", header= TRUE)
df4 <- read.csv2("2009.07.04.csv", sep = ",", header= TRUE)
df5 <- read.csv2("2009.07.05.csv", sep = ",", header= TRUE)
# The above CSV files are zipped and uploaded in the following Google Drive link:
# https://drive.google.com/drive/folders/15vRcoLbTE374qC9pTLfsZtF614zIRXrG?usp=sharing

# Load package
library(igraph)

# Make the network
net1 <- graph_from_data_frame(d=df1, directed=TRUE)
print(net1, e=TRUE, v=TRUE)
net2 <- graph_from_data_frame(d=df2, directed=TRUE)
print(net2, e=TRUE, v=TRUE)
net3 <- graph_from_data_frame(d=df3, directed=TRUE)
print(net3, e=TRUE, v=TRUE)
net4 <- graph_from_data_frame(d=df4, directed=TRUE)
print(net4, e=TRUE, v=TRUE)
net5 <- graph_from_data_frame(d=df5, directed=TRUE)
print(net5, e=TRUE, v=TRUE)

#################### Part 2 ####################

# Number of vertices
verticesFrame <- data.frame("VerticesNumber"=c(vcount(net1), vcount(net2), vcount(net3), vcount(net4), vcount(net5)), 
                            "When" = c("1st day","2nd day", "3rd day", "4th day", "5th day"))

# Number of edges
edgesFrame <- data.frame("EdgesNumber"=c(ecount(net1), ecount(net2), ecount(net3), ecount(net4), ecount(net5)), 
                            "When" = c("1st day","2nd day", "3rd day", "4th day", "5th day"))

# Diameter
diameterFrame <- data.frame("Diameter"=c(diameter(net1), diameter(net2), diameter(net3), diameter(net4), diameter(net5)), 
                         "When" = c("1st day","2nd day", "3rd day", "4th day", "5th day"))

# Average in-degree
in_deg1 <- mean(degree(net1, mode="in"))
in_deg2 <- mean(degree(net2, mode="in"))
in_deg3 <- mean(degree(net3, mode="in"))
in_deg4 <- mean(degree(net4, mode="in"))
in_deg5 <- mean(degree(net5, mode="in"))
inDegFrame <- data.frame("InDegree"=c(in_deg1, in_deg2, in_deg3, in_deg4, in_deg5), 
                            "When" = c("1st day","2nd day", "3rd day", "4th day", "5th day"))

# Average out-degree
out_deg1 <- mean(degree(net1, mode="out"))
out_deg2 <- mean(degree(net2, mode="out"))
out_deg3 <- mean(degree(net3, mode="out"))
out_deg4 <- mean(degree(net4, mode="out"))
out_deg5 <- mean(degree(net5, mode="out"))
outDegFrame <- data.frame("OutDegree"=c(out_deg1, out_deg2, out_deg3, out_deg4, out_deg5), 
                          "When" = c("1st day","2nd day", "3rd day", "4th day", "5th day"))

# Select the colors that will be used
library(RColorBrewer)
# All palette available from RColorBrewer
display.brewer.all()
# Select the first 5 colors in the Set1 palette
cols <- brewer.pal(n=5, name="Set1")
# Colours contain the names of fivr different colors
# Create a color vector corresponding to levels in the When variable in the created data frames
cols_t1 <- cols[verticesFrame$When]
cols_t2 <- cols[edgesFrame$When]
cols_t3 <- cols[diameterFrame$When]
cols_t4 <- cols[inDegFrame$When]
cols_t5 <- cols[outDegFrame$When]

# Vertices plot
plot(verticesFrame$VerticesNumber,  type="b", xlab="Days", ylab="Vertices", col=cols_t1, pch=16, main="Number of network vertices per day")
legend("topright",legend=c("1st of July", "2nd of July", "3rd of July", "4th of July", "5th of July"),
       col=rep(cols,times=2),pch=rep(c(16,18),each=4),bty="n", ncol=2,cex=0.7,pt.cex=0.7)

# Edges plot
plot(edgesFrame$EdgesNumber,  type="b", xlab="Days", ylab="Edges", col=cols_t2, pch=16, main="Number of network edges per day")
legend("topright",legend=c("1st of July", "2nd of July", "3rd of July", "4th of July", "5th of July"),
       col=rep(cols,times=2),pch=rep(c(16,18),each=4),bty="n",ncol=2,cex=0.7,pt.cex=0.7)

# Diameter plot
plot(diameterFrame$Diameter,  type="b", xlab="Days", ylab="Diameter", col=cols_t3, pch=16, main="Network diameter per day")
legend("top",legend=c("1st of July", "2nd of July", "3rd of July", "4th of July", "5th of July"),
       col=rep(cols,times=2),pch=rep(c(16,18),each=4),bty="n",ncol=2,cex=0.7,pt.cex=0.7)

# Average in-degree plot
plot(inDegFrame$InDegree,  type="b", xlab="Days", ylab="In-degree level", col=cols_t4, pch=16, main="Average network in-degree per day")
legend("topleft",legend=c("1st of July", "2nd of July", "3rd of July", "4th of July", "5th of July"),
       col=rep(cols,times=2),pch=rep(c(16,18),each=4),bty="n",ncol=2,cex=0.7,pt.cex=0.7)

# Average out-degree plot
plot(outDegFrame$OutDegree,  type="b", xlab="Days", ylab="Out-degree level", col=cols_t5, pch=16, main="Average network out-degree per day")
legend("topleft",legend=c("1st of July", "2nd of July", "3rd of July", "4th of July", "5th of July"),
       col=rep(cols,times=2),pch=rep(c(16,18),each=4),bty="n",ncol=2,cex=0.7,pt.cex=0.7)

#################### Part 3 ####################

# Top 10 Twitter users (in-degree) 
inDegr1 <- sort(degree(net1, mode="in"), decreasing = TRUE)
top_10_inDegr1 <- head(inDegr1, n=10)
print(top_10_inDegr1)
inDegr2 <- sort(degree(net2, mode="in"), decreasing = TRUE)
top_10_inDegr2 <- head(inDegr2, n=10)
print(top_10_inDegr2)
inDegr3 <- sort(degree(net3, mode="in"), decreasing = TRUE)
top_10_inDegr3 <- head(inDegr3, n=10)
print(top_10_inDegr3)
inDegr4 <- sort(degree(net4, mode="in"), decreasing = TRUE)
top_10_inDegr4 <- head(inDegr4, n=10)
print(top_10_inDegr4)
inDegr5 <- sort(degree(net5, mode="in"), decreasing = TRUE)
top_10_inDegr5 <- head(inDegr5, n=10)
print(top_10_inDegr5)

# Top 10 Twitter users (out-degree) 
outDegr1 <- sort(degree(net1, mode="out"), decreasing = TRUE)
top_10_outDegr1 <- head(outDegr1, n=10)
print(top_10_outDegr1)
outDegr2 <- sort(degree(net2, mode="out"), decreasing = TRUE)
top_10_outDegr2 <- head(outDegr2, n=10)
print(top_10_outDegr2)
outDegr3 <- sort(degree(net3, mode="out"), decreasing = TRUE)
top_10_outDegr3 <- head(outDegr3, n=10)
print(top_10_outDegr3)
outDegr4 <- sort(degree(net4, mode="out"), decreasing = TRUE)
top_10_outDegr4 <- head(outDegr4, n=10)
print(top_10_outDegr4)
outDegr5 <- sort(degree(net5, mode="out"), decreasing = TRUE)
top_10_outDegr5 <- head(outDegr5, n=10)
print(top_10_outDegr5)

# Top 10 Twitter users (PageRank) 

library(magrittr)
pr1 <- net1 %>%
  page.rank(directed = TRUE) %>%
  use_series("vector") %>%
  sort(decreasing = TRUE) %>%
  as.matrix %>%
  set_colnames("page.rank")
top_10_PageRank1 <- head(pr1, n=10)
print(top_10_PageRank1)

pr2 <- net2 %>%
  page.rank(directed = TRUE) %>%
  use_series("vector") %>%
  sort(decreasing = TRUE) %>%
  as.matrix %>%
  set_colnames("page.rank")
top_10_PageRank2 <- head(pr2, n=10)
print(top_10_PageRank2)

pr3 <- net3 %>%
  page.rank(directed = TRUE) %>%
  use_series("vector") %>%
  sort(decreasing = TRUE) %>%
  as.matrix %>%
  set_colnames("page.rank")
top_10_PageRank3 <- head(pr3, n=10)
print(top_10_PageRank3)


pr4 <- net4 %>%
  page.rank(directed = TRUE) %>%
  use_series("vector") %>%
  sort(decreasing = TRUE) %>%
  as.matrix %>%
  set_colnames("page.rank")
top_10_PageRank4 <- head(pr4, n=10)
print(top_10_PageRank4)

pr5 <- net5 %>%
  page.rank(directed = TRUE) %>%
  use_series("vector") %>%
  sort(decreasing = TRUE) %>%
  as.matrix %>%
  set_colnames("page.rank")
top_10_PageRank5 <- head(pr5, n=10)
print(top_10_PageRank5)

#################### Part 4 ####################

######  Question 1: Algorithms ###### 
      
      # Make network graph undirected
      net1_undir <- as.undirected(net1)
      net2_undir <- as.undirected(net2)
      net3_undir <- as.undirected(net3)
      net4_undir <- as.undirected(net4)
      net5_undir <- as.undirected(net5)
      
      # Find communities with fast greedy clustering
      communities_fast_greedy1 <- cluster_fast_greedy(net1_undir)
      communities_fast_greedy2 <- cluster_fast_greedy(net2_undir)
      communities_fast_greedy3 <- cluster_fast_greedy(net3_undir)
      communities_fast_greedy4 <- cluster_fast_greedy(net4_undir)
      communities_fast_greedy5 <- cluster_fast_greedy(net5_undir)
      
      # Find communities with infomap clustering
      communities_infomap1 <- cluster_infomap(net1_undir)
      communities_infomap2 <- cluster_infomap(net2_undir)
      communities_infomap3 <- cluster_infomap(net3_undir)
      communities_infomap4 <- cluster_infomap(net4_undir)
      communities_infomap5 <- cluster_infomap(net5_undir)
      
      # Find communities with louvain clustering
      communities_louvain1 <- cluster_louvain(net1_undir)
      communities_louvain2 <- cluster_louvain(net2_undir)
      communities_louvain3 <- cluster_louvain(net3_undir)
      communities_louvain4 <- cluster_louvain(net4_undir)
      communities_louvain5 <- cluster_louvain(net5_undir)
      
      # Compare fast greedy communities with louvain clustering
      compare(communities_fast_greedy1, communities_louvain1)
      compare(communities_fast_greedy2, communities_louvain2)
      compare(communities_fast_greedy3, communities_louvain3)
      compare(communities_fast_greedy4, communities_louvain4)
      compare(communities_fast_greedy5, communities_louvain5)

###### Question 2: Pick a random user and examine the evolvement of the communities this user belongs to ###### 
      
      # Select a user who appears in all 5 graphs (e.g. "egomonics")
      commonUsers<-Reduce(intersect, list(communities_louvain1$names, communities_louvain2$names, 
                             communities_louvain3$names, communities_louvain4$names,
                             communities_louvain5$names))
      
      # Search the community ID of the communitites in which the user "egomonics" belongs to
      object1 <- membership(communities_louvain1)
      object2 <- membership(communities_louvain2)
      object3 <- membership(communities_louvain3)
      object4 <- membership(communities_louvain4)
      object5 <- membership(communities_louvain5)
      
      library(reshape2)
      members1 <- melt(as.matrix(object1), varnames=c("TwitterUser","Number"),value.name="CommunityID")
      members1$Number<-NULL
      members2 <- melt(as.matrix(object2), varnames=c("TwitterUser","Number"),value.name="CommunityID")
      members2$Number<-NULL
      members3 <- melt(as.matrix(object3), varnames=c("TwitterUser","Number"),value.name="CommunityID")
      members3$Number<-NULL
      members4 <- melt(as.matrix(object4), varnames=c("TwitterUser","Number"),value.name="CommunityID")
      members4$Number<-NULL
      members5 <- melt(as.matrix(object5), varnames=c("TwitterUser","Number"),value.name="CommunityID")
      members5$Number<-NULL
      
      userCommunityID1<-subset(members1, TwitterUser==commonUsers[2])
      userCommunityID2<-subset(members2, TwitterUser==commonUsers[2])
      userCommunityID3<-subset(members3, TwitterUser==commonUsers[2])
      userCommunityID4<-subset(members4, TwitterUser==commonUsers[2])
      userCommunityID5<-subset(members5, TwitterUser==commonUsers[2])
      
      # 1st criterio for examination: Number of vertices among communities

            communityVertices <- data.frame("VerticesNumber"=c(length(which(communities_louvain1$membership==userCommunityID1$CommunityID)), 
                                                                      length(which(communities_louvain2$membership==userCommunityID2$CommunityID)), 
                                                                      length(which(communities_louvain3$membership==userCommunityID3$CommunityID)), 
                                                                      length(which(communities_louvain4$membership==userCommunityID4$CommunityID)),
                                                                      length(which(communities_louvain5$membership==userCommunityID5$CommunityID))), 
                                             "When" = c("1st day","2nd day", "3rd day", "4th day", "5th day"))
            # Select the colors that will be used
            library(RColorBrewer)
            # All palette available from RColorBrewer
            display.brewer.all()
            # Select the first 5 colors in the Set1 palette
            cols2 <- brewer.pal(n=5, name="Set1")
            # Colours contain the names of fivr different colors
            # Create a color vector corresponding to levels in the When variable in the created data frames
            cols_t2 <- cols2[communityVertices$When]
            
            # Vertices plot for the communities of user egomonics
            plot(communityVertices$VerticesNumber,  type="b", xlab="Days", ylab="Vertices", col=cols_t2, pch=16, main="Number of vertices per community")
            legend("top",legend=c("1st community", "2nd community", "3rd community", "4th community", "5th community"),
                   col=rep(cols2,times=2),pch=rep(c(16,18),each=4),bty="n", ncol=2,cex=0.7,pt.cex=0.7)

      # 2nd criterio for examination: Similarity of users among communities
      
            community1Users <- subset(members1, CommunityID == 76644) 
            community1Users$TwitterUser<- as.character(community1Users$TwitterUser)
            community2Users <- subset(members2, CommunityID == 56224) 
            community2Users$TwitterUser<- as.character(community2Users$TwitterUser)
            community3Users <- subset(members3, CommunityID == 28998) 
            community3Users$TwitterUser<- as.character(community3Users$TwitterUser)
            community4Users <- subset(members4, CommunityID == 22543) 
            community4Users$TwitterUser<- as.character(community4Users$TwitterUser)
            community5Users <- subset(members5, CommunityID == 16117) 
            community5Users$TwitterUser<- as.character(community5Users$TwitterUser)
            
            # Count the number of users who simultaneously belonged to 1,2,3,4 or 5 communities
            
            combinedDf <- data.frame(table(c(community1Users$TwitterUser, community2Users$TwitterUser, community3Users$TwitterUser, community4Users$TwitterUser, community5Users$TwitterUser)))
            names(combinedDf) <- c("TwitterUser", "Matches")
            
            library(dplyr)
            usersParticipatingIn1Communities <- nrow(as.data.frame (combinedDf %>% group_by(Matches) %>% filter(Matches==1)))
            usersParticipatingIn2Communities <- nrow(as.data.frame (combinedDf %>% group_by(Matches) %>% filter(Matches==2)))
            usersParticipatingIn3Communities <- nrow(as.data.frame (combinedDf %>% group_by(Matches) %>% filter(Matches==3)))
            usersParticipatingIn4Communities <- nrow(as.data.frame (combinedDf %>% group_by(Matches) %>% filter(Matches==4)))
            usersParticipatingIn5Communities <- nrow(as.data.frame (combinedDf %>% group_by(Matches) %>% filter(Matches==5)))
            membersPercentageSimilarity <- data.frame("Percentage"=c((usersParticipatingIn1Communities/nrow(combinedDf))*100,
                                                               (usersParticipatingIn2Communities/nrow(combinedDf))*100,
                                                               (usersParticipatingIn3Communities/nrow(combinedDf))*100,
                                                               (usersParticipatingIn4Communities/nrow(combinedDf))*100,
                                                               (usersParticipatingIn5Communities/nrow(combinedDf))*100), 
                                                              "WhichCommunity" = c("1st","2nd", "3rd", "4th", "5th"))
            # Plot the results
              # Select the colors that will be used
              library(RColorBrewer)
              # All palette available from RColorBrewer
              display.brewer.all()
              # Select the first 5 colors in the Set1 palette
              cols3 <- brewer.pal(n=5, name="Set1")
              # Create a color vector corresponding to levels in the When variable in the created data frames
              cols_t3 <- cols3[membersPercentageSimilarity$WhichCommunity]
              # Plot and depict if the same users existed in the communities in which user egomonics belongs to
              plot(membersPercentageSimilarity$Percentage,  type="b", xlab="Users involvement in each community", ylab="Percentage", col=cols_t2, pch=16, main="Users' similarity among communities")
              legend("topright",legend=c("1st community", "2nd community", "3rd community", "4th community", "5th community"),
                     col=rep(cols3,times=2),pch=rep(c(16,18),each=4),bty="n", ncol=2,cex=0.7,pt.cex=0.7)

###### Question 3: Visualize in a graph the communities with different colours ###### 
            
# Communities of undirected network 1
    # Get the sizes of each community
    community_size1 <- sizes(communities_louvain1)
    # Some mid-size communities
    in_mid_community1 <- unlist(communities_louvain1[community_size1 > 50 & community_size1 < 90])
    # Induce a subgraph of graph using in_mid_community
    subgraph1_directed <- induced.subgraph(net1, in_mid_community1)
    subgraph1_undirected<- as.undirected(subgraph1_directed)
    subgraph_louvain1 <- cluster_louvain(subgraph1_undirected)
    # Does the edge cross betwen commmunities?
    is_crossing1 <- crossing(net1, communities = communities_louvain1)
    # Set edge linetype: solid for crossings, dotted otherwise 
     E(net1)$lty <- ifelse(is_crossing1, "solid", "dotted")
    # Plot those mid-size communities
    plot(subgraph1_directed, vertex.color=rainbow(16, alpha=1)[subgraph_louvain1$membership], vertex.label=NA, edge.arrow.size=.2, vertex.size=10,
         margin = 0,  coords = layout_with_fr(subgraph1_directed), edge.arrow.width = 0.8, edge.arrow.size = 0.2, lty=E(net1)$lty, main="Communities of 1st day")
    
# Communities of undirected network 2
    # Get the sizes of each community
    community_size2 <- sizes(communities_louvain2)
    # Some mid-size communities
    in_mid_community2 <- unlist(communities_louvain2[community_size2 > 50 & community_size2 < 90])
    # Induce a subgraph of graph using in_mid_community
    subgraph2_directed <- induced.subgraph(net2, in_mid_community2)
    subgraph2_undirected<- as.undirected(subgraph2_directed)
    subgraph_louvain2 <- cluster_louvain(subgraph2_undirected)
    # Does the edge cross betwen commmunities?
    is_crossing2 <- crossing(net2, communities = communities_louvain2)
    # Set edge linetype: solid for crossings, dotted otherwise 
    E(net2)$lty <- ifelse(is_crossing2, "solid", "dotted")
    # Plot those mid-size communities
    plot(subgraph2_directed, vertex.color=rainbow(16, alpha=1)[subgraph_louvain2$membership], vertex.label=NA, edge.arrow.size=.2, vertex.size=10,
         margin = 0,  coords = layout_with_fr(subgraph2_directed), edge.arrow.width = 0.8, edge.arrow.size = 0.2, lty= E(net2)$lty, main="Communities of 2nd day")
    
# Communities of undirected network 3
    # Get the sizes of each community
    community_size3 <- sizes(communities_louvain3)
    # Some mid-size communities
    in_mid_community3 <- unlist(communities_louvain3[community_size3 > 50 & community_size3 < 90])
    # Induce a subgraph of graph using in_mid_community
    subgraph3_directed <- induced.subgraph(net3, in_mid_community3)
    subgraph3_undirected<- as.undirected(subgraph3_directed)
    subgraph_louvain3 <- cluster_louvain(subgraph3_undirected)
    # Does the edge cross betwen commmunities?
    is_crossing3 <- crossing(net3, communities = communities_louvain3)
    # Set edge linetype: solid for crossings, dotted otherwise 
    E(net3)$lty <- ifelse(is_crossing3, "solid", "dotted")
    # Plot those mid-size communities
    plot(subgraph3_directed, vertex.color=rainbow(16, alpha=1)[subgraph_louvain3$membership], vertex.label=NA, edge.arrow.size=.2, vertex.size=10,
         margin = 0,  coords = layout_with_fr(subgraph3_directed), edge.arrow.width = 0.8, edge.arrow.size = 0.2, lty=E(net3)$lty, main="Communities of 3rd day")
    
# Communities of undirected network 4
    # Get the sizes of each community
    community_size4 <- sizes(communities_louvain4)
    # Some mid-size communities
    in_mid_community4 <- unlist(communities_louvain4[community_size4 > 50 & community_size4 < 90])
    # Induce a subgraph of graph using in_mid_community
    subgraph4_directed <- induced.subgraph(net4, in_mid_community4)
    subgraph4_undirected<- as.undirected(subgraph4_directed)
    subgraph_louvain4 <- cluster_louvain(subgraph4_undirected)
    # Does the edge cross betwen commmunities?
    is_crossing4 <- crossing(net4, communities = communities_louvain4)
    # Set edge linetype: solid for crossings, dotted otherwise 
    E(net4)$lty <- ifelse(is_crossing4, "solid", "dotted")
    # Plot those mid-size communities
    plot(subgraph4_directed, vertex.color=rainbow(16, alpha=1)[subgraph_louvain4$membership], vertex.label=NA, edge.arrow.size=.2, vertex.size=10,
         margin = 0,  coords = layout_with_fr(subgraph4_directed), edge.arrow.width = 0.8, edge.arrow.size = 0.2, lty=E(net4)$lty, main="Communities of 4th day")
    
# Communities of undirected network 5
    # Get the sizes of each community
    community_size5 <- sizes(communities_louvain5)
    # Some mid-size communities
    in_mid_community5 <- unlist(communities_louvain5[community_size5 > 50 & community_size5 < 90])
    # Induce a subgraph of graph using in_mid_community
    subgraph5_directed <- induced.subgraph(net5, in_mid_community5)
    subgraph5_undirected<- as.undirected(subgraph5_directed)
    subgraph_louvain5 <- cluster_louvain(subgraph5_undirected)
    # Does the edge cross betwen commmunities?
    is_crossing5 <- crossing(net5, communities = communities_louvain5)
    # Set edge linetype: solid for crossings, dotted otherwise 
    E(net5)$lty <- ifelse(is_crossing5, "solid", "dotted")
    # Plot those mid-size communities
    plot(subgraph5_directed, vertex.color=rainbow(16, alpha=1)[subgraph_louvain5$membership], vertex.label=NA, edge.arrow.size=.2, vertex.size=10,
         margin = 0,  coords = layout_with_fr(subgraph5_directed), edge.arrow.width = 0.8, edge.arrow.size = 0.2, lty=E(net5)$lty, main="Communities of 5th day")
