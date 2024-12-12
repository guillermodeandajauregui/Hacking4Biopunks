# excercise 2: MI filtering and network construction----

library(tidyverse)
library(igraph)
library(tidygraph)
library(ggraph)

# read data ----

mi_matrix <- read_delim(file = "mi_matrix.txt", delim = "\t")

y <- 
  mi_matrix |> 
  select(-1) |> 
  as.matrix()

rownames(y) <- mi_matrix$gene  

# diagonal by definition should be zero (avoid self loops) ----
diag(y) <- 0

################################################################################
# binarize  using a threshold ----
################################################################################

threshold <-quantile(x = y, probs = 0.99)

z <- ifelse(y >= threshold, 1, 0)

# and make a network out of it ----

g <- graph_from_adjacency_matrix(z, mode = "undirected")

mi_layout = layout_with_kk(g)

plot(g, layout = mi_layout, vertex.size = 1, vertex.label="")

#write_graph(g, file = "example_graph.graphml", format = "graphml")

# what if we try other thresholds?


################################################################################
# binarize  using a local threshold ----

# time intensive!!!
################################################################################


# let's define the disparity function ----

# Function to calculate disparity
calculate_disparity <- function(g, weight_attr = "weight") {
  # Create a new graph to store the results
  B <- g
  
  # Get the number of vertices
  vertices <- V(g)
  
  for (u in vertices) {
    # Get neighbors of node u
    neighbors_u <- neighbors(g, u)
    k <- length(neighbors_u)

        # Continue only if node has more than one neighbor
    if (k > 1) {
      # Sum of weights of edges connected to node u
      # Get the edges incident to node u
      edges_u <- E(g)[.inc(u)]
      # Sum of weights of edges connected to node u
      
      weights <- edge_attr(graph = g, name = weight_attr, index = edges_u)

      sum_w <- sum(weights)  # Access edge weights directly

            # Iterate through each neighbor v of u
      for (v in neighbors_u) {
        # Get the edge index between u and v
        edge_uv <- E(g, P = c(u, v))
        
        # Get the weight of the edge between u and v
        w <- edge_attr(graph = g, name = weight_attr, index = edge_uv)
        # Calculate p_ij
        p_ij <- w / sum_w
        
        # Calculate alpha_ij using integration
        alpha_ij <- 1 - (k - 1) * integrate(function(x) (1 - x)^(k - 2), 0, p_ij)$value
        
        # Set alpha value for the edge
        E(B)[edge_uv]$alpha <- round(alpha_ij, 4)
      }
    }
  }
  
  return(B)
}

set.seed(725)
{g_er <- erdos.renyi.game(n = 50, p.or.m = 0.5, type = "gnp", directed = "F")%>%
  set_edge_attr("weight", value = runif(ecount(.), 1, 10))
}
# Calculate disparity
g_er <- calculate_disparity(g_er)

mi_layout = layout_with_kk(g_er)
plot(g_er, layout = mi_layout, vertex.size = 1, vertex.label="")

g_er_filtered <- delete_edges(graph = g_er, edges = E(g_er)[alpha > 0.4])
plot(g_er_filtered, layout = mi_layout, vertex.size = 1, vertex.label="")


################################################################################
# binarize  using a dpi ----

# time intensive!!!
################################################################################

# Function to apply Data Processing Inequality filter with a threshold
dpi_filter <- function(mi_matrix, threshold = 0.15) {
  n <- nrow(mi_matrix)
  
  # Loop over all triplets (i, j, k)
  for (i in 1:(n-2)) {
    for (j in (i+1):(n-1)) {
      for (k in (j+1):n) {
        mi_ij <- mi_matrix[i, j]
        mi_ik <- mi_matrix[i, k]
        mi_jk <- mi_matrix[j, k]
        
        # Remove edges based on DPI and threshold
        if (mi_ik < threshold || mi_ik > min(mi_ij, mi_jk)) {
          mi_matrix[i, k] <- 0  # Remove indirect edge i-k
          mi_matrix[k, i] <- 0
        }
        if (mi_jk < threshold || mi_jk > min(mi_ij, mi_ik)) {
          mi_matrix[j, k] <- 0  # Remove indirect edge j-k
          mi_matrix[k, j] <- 0
        }
      }
    }
  }
  
  return(mi_matrix)
}

m_er <- as_adjacency_matrix(g_er, attr = "weight")

filtered_mi_matrix <- dpi_filter(m_er, threshold = 0.15)


mi_layout = layout_with_kk(g_er)
plot(g_er, layout = mi_layout, vertex.size = 1, vertex.label="")

g_er_filtered <- graph_from_adjacency_matrix(filtered_mi_matrix, mode = "undirected", weighted = T)
plot(g_er_filtered, layout = mi_layout, vertex.size = 1, vertex.label="")
