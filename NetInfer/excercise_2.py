import pandas as pd
import numpy as np
from sklearn.metrics import mutual_info_score
from sklearn.preprocessing import KBinsDiscretizer
import seaborn as sns
import matplotlib.pyplot as plt

# read data ----
os.chdir('~/GITS/workshop4ukraine/')  # Replace with your actual path

x = pd.read_csv('gene_data.txt', sep='\t')

# make it a matrix ----
y = x.iloc[:, 2:].values  # select all columns except the first two
rownames = x['hgnc_symbol'].values  # store the rownames

# discretize data ----
# Function to discretize each row
def discretize_row(row, n_bins):
    # Use 'quantile' strategy for equal frequency discretization
    discretizer = KBinsDiscretizer(n_bins=n_bins, encode='ordinal', strategy='quantile')
    row_discretized = discretizer.fit_transform(row.reshape(-1, 1)).flatten()  # Discretize
    return row_discretized

# Apply discretization to each row, adjust bin size
n_bins = round(len(y[0]) ** (1/3))  # cube root rule for the number of bins
y_discretized = [discretize_row(row, n_bins=n_bins) for row in y]

# calculate MI ----
mi_matrix = np.zeros((len(y_discretized), len(y_discretized)))

for i in range(len(y_discretized)):
    for j in range(len(y_discretized)):
        mi_matrix[i, j] = mutual_info_score(y_discretized[i], y_discretized[j])

# plot MI matrix ----
sns.heatmap(mi_matrix, xticklabels=rownames, yticklabels=rownames, cmap="viridis")
plt.title('Mutual Information Matrix')
plt.show()
