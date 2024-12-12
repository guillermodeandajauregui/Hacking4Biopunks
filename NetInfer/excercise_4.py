import networkx as nx
from networkx.algorithms import community
import matplotlib.pyplot as plt
import numpy as np
from pygraphml import GraphMLParser
import pandas as pd

################################################################################
# read networks ----
################################################################################

# read a graphml file ----
parser = GraphMLParser()
g = nx.read_graphml('example_graph.graphml')

################################################################################
# analyze paths ----
################################################################################

# Shortest path matrix
path_matrix = dict(nx.all_pairs_shortest_path_length(g))
# Convert to a DataFrame for better viewing
path_df = pd.DataFrame(path_matrix)
print(path_df.head())

# Average shortest path length
avg_shortest_path = nx.average_shortest_path_length(g)
print("Average shortest path length:", avg_shortest_path)

# Diameter
my_diameter = nx.diameter(g)
print("Diameter:", my_diameter)

################################################################################
# connected components ----
################################################################################

components_g = list(nx.connected_components(g))
num_components = nx.number_connected_components(g)
component_sizes = [len(c) for c in components_g]
print("Number of components:", num_components)
print("Component sizes:", component_sizes)

# Degree ----
degree_dict = dict(g.degree())
nx.set_node_attributes(g, degree_dict, 'degree')
print("Degree of nodes:", degree_dict)

# If the graph is directed, calculate in-degree and out-degree
if nx.is_directed(g):
    in_degree_dict = dict(g.in_degree())
    out_degree_dict = dict(g.out_degree())
    nx.set_node_attributes(g, in_degree_dict, 'in_degree')
    nx.set_node_attributes(g, out_degree_dict, 'out_degree')

################################################################################
# strength (weighted degree) ----
################################################################################

# Strength (sum of weights of edges incident to each node)
strength_dict = dict(g.degree(weight='weight'))
nx.set_node_attributes(g, strength_dict, 'strength')
print("Strength of nodes:", strength_dict)

################################################################################
# betweenness centrality ----
################################################################################

# Node betweenness centrality
betweenness_centrality = nx.betweenness_centrality(g)
nx.set_node_attributes(g, betweenness_centrality, 'betweenness_centrality')

# Edge betweenness centrality
edge_betweenness = nx.edge_betweenness_centrality(g)
nx.set_edge_attributes(g, edge_betweenness, 'edge_betweenness')

################################################################################
# clustering coefficient ----
################################################################################

# Local clustering coefficient
clustering_coeffs = nx.clustering(g)
nx.set_node_attributes(g, clustering_coeffs, 'clustering_coefficient')

# Global clustering coefficient
cc_global = nx.transitivity(g)
print("Global clustering coefficient:", cc_global)

################################################################################
# community detection ----
################################################################################

# Louvain method for community detection
import community as community_louvain
partition_louvain = community_louvain.best_partition(g)
nx.set_node_attributes(g, partition_louvain, 'community_louvain')

# Other community detection algorithms
partition_infomap = list(nx.algorithms.community.label_propagation.label_propagation_communities(g))
partition_walktrap = list(nx.algorithms.community.asyn_fluid.asyn_fluidc(g, k=5)) # Example
partition_fstgred = community.greedy_modularity_communities(g)

# Assign the detected communities to the nodes
for idx, community in enumerate(partition_walktrap):
    for node in community:
        g.nodes[node]['community_walktrap'] = idx

################################################################################
# plot communities ----
################################################################################

# Visualize the network and communities
pos = nx.spring_layout(g)  # Fruchterman-Reingold layout

# Plot Louvain community detection
plt.figure(figsize=(10, 10))
nx.draw(g, pos, node_color=list(partition_louvain.values()), cmap=plt.get_cmap('viridis'), with_labels=False, node_size=50)
plt.title("Louvain Community Detection")
plt.show()

# Check if an edge is inter- or intra-community for Louvain
crossing_louvain = nx.edge_boundary(g, partition_louvain)
nx.set_edge_attributes(g, {edge: True for edge in crossing_louvain}, 'crossing_louvain')

################################################################################
# Visualize using different layouts and libraries ----
################################################################################

# Plot using the Kamada-Kawai layout
plt.figure(figsize=(10, 10))
nx.draw(g, pos, node_color=list(partition_louvain.values()), cmap=plt.get_cmap('viridis'), with_labels=False, node_size=50, edge_color='grey')
plt.title("Kamada-Kawai Layout with Louvain Community")
plt.show()

# Plot using circular layout
pos_circle = nx.circular_layout(g)
plt.figure(figsize=(10, 10))
nx.draw(g, pos_circle, node_color=list(partition_louvain.values()), cmap=plt.get_cmap('viridis'), with_labels=False, node_size=50, edge_color='grey')
plt.title("Circular Layout with Louvain Community")
plt.show()

# Plot using Fruchterman-Reingold layout (default spring layout)
plt.figure(figsize=(10, 10))
nx.draw(g, pos, node_color=list(partition_louvain.values()), cmap=plt.get_cmap('viridis'), with_labels=False, node_size=50, edge_color='grey')
plt.title("Fruchterman-Reingold Layout with Louvain Community")
plt.show()
