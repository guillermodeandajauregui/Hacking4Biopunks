import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from scipy.integrate import quad

# read data ----
mi_matrix = pd.read_csv('mi_matrix.txt', sep='\t')

# Convert to numpy matrix and set row names
y = mi_matrix.iloc[:, 1:].values
genes = mi_matrix['gene'].values

# Diagonal by definition should be zero (avoid self-loops)
np.fill_diagonal(y, 0)

################################################################################
# binarize using a threshold ----
################################################################################

threshold = np.quantile(y, 0.9)
z = np.where(y >= threshold, 1, 0)

# make a network out of it ----
G = nx.from_numpy_matrix(z, create_using=nx.Graph)
pos = nx.kamada_kawai_layout(G)
nx.draw(G, pos, node_size=50, with_labels=False)
plt.show()

################################################################################
# binarize using a local threshold (disparity function) ----

# Function to calculate disparity
def calculate_disparity(G, weight_attr='weight'):
    B = G.copy()
    
    for u in G.nodes():
        neighbors_u = list(G.neighbors(u))
        k = len(neighbors_u)
        
        if k > 1:
            # Sum of weights of edges connected to node u
            sum_w = sum(abs(G[u][v][weight_attr]) for v in neighbors_u)
            
            # Iterate through each neighbor v of u
            for v in neighbors_u:
                w = abs(G[u][v][weight_attr])
                p_ij = float(w) / sum_w
                
                # Calculate alpha_ij using integration
                alpha_ij = 1 - (k - 1) * quad(lambda x: (1 - x)**(k - 2), 0, p_ij)[0]
                B[u][v]['alpha'] = round(alpha_ij, 4)
    
    return B

# Example random graph with weights
np.random.seed(725)
G_er = nx.erdos_renyi_graph(50, 0.5)
for (u, v) in G_er.edges():
    G_er[u][v]['weight'] = np.random.uniform(1, 10)

# Calculate disparity
G_er = calculate_disparity(G_er)

# Plot original network
pos = nx.kamada_kawai_layout(G_er)
nx.draw(G_er, pos, node_size=50, with_labels=False)
plt.show()

# Filter edges with alpha > 0.4
edges_to_remove = [(u, v) for u, v, alpha in G_er.edges(data='alpha') if alpha > 0.4]
G_er_filtered = G_er.copy()
G_er_filtered.remove_edges_from(edges_to_remove)

# Plot filtered network
nx.draw(G_er_filtered, pos, node_size=50, with_labels=False)
plt.show()

################################################################################
# binarize using a dpi (Data Processing Inequality) filter ----

# Function to apply DPI filter with a threshold
def dpi_filter(mi_matrix, threshold=0.15):
    n = mi_matrix.shape[0]
    
    for i in range(n - 2):
        for j in range(i + 1, n - 1):
            for k in range(j + 1, n):
                mi_ij = mi_matrix[i, j]
                mi_ik = mi_matrix[i, k]
                mi_jk = mi_matrix[j, k]
                
                # Remove edges based on DPI and threshold
                if mi_ik < threshold or mi_ik > min(mi_ij, mi_jk):
                    mi_matrix[i, k] = 0
                    mi_matrix[k, i] = 0
                if mi_jk < threshold or mi_jk > min(mi_ij, mi_ik):
                    mi_matrix[j, k] = 0
                    mi_matrix[k, j] = 0
    
    return mi_matrix

# Convert G_er to adjacency matrix with weights
m_er = nx.to_numpy_matrix(G_er, weight='weight')

# Apply DPI filter
filtered_mi_matrix = dpi_filter(m_er, threshold=0.15)

# Create a new graph from the filtered matrix
G_er_filtered = nx.from_numpy_matrix(filtered_mi_matrix, create_using=nx.Graph)

# Plot the filtered network
nx.draw(G_er_filtered, pos, node_size=50, with_labels=False)
plt.show()
