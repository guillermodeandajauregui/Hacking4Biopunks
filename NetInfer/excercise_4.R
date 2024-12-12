#excercise 4: network analysis ----

library(igraph)

################################################################################
#read networks ----
################################################################################


#read a graphml ----

g <- read_graph(file = "example_graph.graphml", "graphml")




################################################################################
#analyze paths ----
################################################################################


E(g)$weight_a <- E(g)$weight
E(g)$weight   <- 1/E(g)$weight

path_matrix <- shortest.paths(g)
path_matrix[1:5, 1:5]

avg_shortest_path <- average.path.length(g)

my_diameter <- diameter(g)

################################################################################
#connected components ----
################################################################################

components_g <- components(g)
components_g$membership
components_g$csize
components_g$no


# degree ----

V(g) 
degree(g)

g <- set.vertex.attribute(graph = g, name = "degree", value = degree(g))
V(g)$degree_2 <- degree(g)

head(get.data.frame(g, what = "vertices"))

V(directed_g)$all_degree <-  degree(directed_g, mode = "all")
V(directed_g)$in_degree <-  degree(directed_g, mode = "in")
V(directed_g)$out_degree <-  degree(directed_g, mode = "out")

head(get.data.frame(directed_g, what = "vertices"))


#strength ---- 


V(g)$strength <- strength(g)
head(get.data.frame(g, what = "vertices"))

# betweenness centrality ----

V(g)$betweenness_centrality <- betweenness(g)
V(g)$betweenness_centrality.estimate <- betweenness.estimate(graph = g, cutoff = 5)
V(g)$betweenness_centrality.estimate.badchoice <- betweenness.estimate(graph = g, cutoff = 3)
head(get.data.frame(g, what = "vertices"))

#edge betweenness ----
E(g)$edge_betweenness <- edge.betweenness(g)
head(get.data.frame(g, what = "edges"))

#clustering coefficient ---- 

V(g)$cc   <- transitivity(graph = g, type = "local", isolates = "zero")
cc_global <- transitivity(graph = g, type = "global")
head(get.data.frame(g, what = "vertices"))
#module detection ---- 

comm.louvain <- cluster_louvain(g)
comm.louvain

#add membership data ---- 
V(g)$comm.louvain <- membership(comm.louvain)
head(get.data.frame(g, what = "vertices"))

comm.infomap <- cluster_infomap(g)
comm.walktrp <- cluster_walktrap(g)
comm.fstgred <- cluster_fast_greedy(g)

#try different algorithms 
V(g)$comm.infomap <- membership(comm.infomap) #random walker based
V(g)$comm.walktrp <- membership(comm.walktrp) #random walker based
V(g)$comm.fstgred <- membership(comm.fstgred) #modularity maximization


head(get.data.frame(g, what = "vertices"))


set.seed(1)
plot(comm.louvain, g, layout = layout_nicely)
set.seed(1)
plot(comm.infomap, g, layout = layout_nicely)
set.seed(1)
plot(comm.walktrp, g, layout = layout_nicely)
set.seed(1)
plot(comm.fstgred, g, layout = layout_nicely)

#is edge inter or intra community? ---- 
E(g)$crossing_louvain <- igraph::crossing(communities = comm.louvain, graph = g)

#some nicer looking plotting 
library(ggraph)

set.seed(1)
ggraph(graph = g, layout = "kk") + 
  geom_edge_link(aes(color = crossing_louvain)) + 
  geom_node_point(aes(color = as.factor(comm.louvain))) + 
  theme_graph() + 
  scale_edge_color_manual("crossing", values = c("grey", "black")) + 
  scale_color_discrete("Community (Louvain)")

set.seed(1)
ggraph(graph = g, layout = "circle") + 
  geom_edge_link(aes(color = crossing_louvain)) + 
  geom_node_point(aes(color = as.factor(comm.louvain))) + 
  theme_graph() + 
  scale_edge_color_manual("crossing", values = c("grey", "black")) + 
  #scale_color_discrete("Community (Louvain)") + 
  scale_color_viridis(name="Community (Louvain)", discrete = T, option = "C")


set.seed(1)
ggraph(graph = g, layout = "fr") + 
  geom_edge_link(aes(color = crossing_louvain)) + 
  geom_node_point(aes(color = as.factor(comm.louvain)), size = 10 ) + 
  theme_graph() + 
  scale_edge_color_manual("crossing", values = c("grey", "black")) + 
  #scale_color_discrete("Community (Louvain)") + 
  scale_color_viridis(name="Community (Louvain)", discrete = T, option = "C")