# -------------------------------------
# 1. Import Libraries
# -------------------------------------
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests

# -------------------------------------
# 2. Load Expression Data and Annotations
# -------------------------------------
# Load expression data
expr_data = pd.read_csv("/datos/rnaseq_demos/GSE67333_norm_counts_FPKM_GRCh38.p13_NCBI.tsv.gz", sep="\t")
expr_data.set_index("GeneID", inplace=True)

# Load sample annotations
annotations = pd.read_csv("/datos/rnaseq_demos/gse67333_labels.txt")
annotations["Condition"] = annotations["Condition"].str[:-1]  # Clean condition labels

# Reorder columns based on annotations
expr_data = expr_data[annotations["SampleID"]]

# -------------------------------------
# 3. Heatmap of Top Variable Genes
# -------------------------------------
# Calculate variance and select top 50 genes
expr_data["Variance"] = expr_data.var(axis=1)
top_genes = expr_data.sort_values("Variance", ascending=False).head(50).drop(columns=["Variance"])

# Plot heatmap
sns.heatmap(top_genes, cmap="viridis", xticklabels=True, yticklabels=False)
plt.title("Heatmap of Top 50 Variable Genes")
plt.show()

# -------------------------------------
# 4. PCA Analysis
# -------------------------------------
# Perform PCA on transposed data (samples as rows)
pca = PCA(n_components=2)
pca_results = pca.fit_transform(expr_data.T)

# Create a dataframe for PCA results
pca_df = pd.DataFrame(pca_results, columns=["PC1", "PC2"], index=expr_data.columns)
pca_df = pca_df.merge(annotations, left_index=True, right_on="SampleID")

# Scatterplot of PCA
sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="Condition", style="Condition", s=100)
plt.title("PCA of Samples")
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)")
plt.legend(bbox_to_anchor=(1, 1))
plt.show()

# -------------------------------------
# 5. Differential Expression Analysis (Simplified Limma-like Approach)
# -------------------------------------
# Define groups
group_ctrl = annotations.loc[annotations["Condition"] == "CTRL", "SampleID"]
group_load = annotations.loc[annotations["Condition"] == "LOAD", "SampleID"]

# Calculate log2 fold change
expr_data["logFC"] = np.log2(expr_data[group_load].mean(axis=1) + 1) - np.log2(expr_data[group_ctrl].mean(axis=1) + 1)

# Perform t-tests (as a proxy for differential expression)
from scipy.stats import ttest_ind

expr_data["P.Value"] = expr_data.apply(
    lambda row: ttest_ind(row[group_load], row[group_ctrl], equal_var=False).pvalue, axis=1
)

# Adjust p-values using Benjamini-Hochberg correction
expr_data["adj.P.Val"] = multipletests(expr_data["P.Value"], method="fdr_bh")[1]

# Save filtered results
filtered = expr_data[(expr_data["logFC"].abs() > 1) & (expr_data["adj.P.Val"] < 0.05)]
filtered.to_csv("/datos/rnaseq_demos/diffex_demo_python.txt", sep="\t")

# -------------------------------------
# 6. Functional Enrichment Analysis (Gene Ontology)
# -------------------------------------
# Load annotation for mapping
annot = pd.read_csv("/datos/rnaseq_demos/Human.GRCh38.p13.annot.tsv.gz", sep="\t")

# Map significant genes to HGNC symbols
diffex_with_hgnc = filtered.reset_index().merge(annot, left_on="GeneID", right_on="GeneID")
gene_list = diffex_with_hgnc["Symbol"].dropna().unique()

# Perform GO enrichment using gseapy
from gseapy import enrichr

# Enrichment analysis
enrich_results = enrichr(gene_list=gene_list, gene_sets="GO_Biological_Process_2021", organism="Human")

# Extract top results
top_go = enrich_results.results.sort_values("Adjusted P-value").head(10)

# Plot top results
sns.barplot(data=top_go, x="-log10(Adjusted P-value)", y="Term")
plt.title("Top GO Enriched Terms")
plt.xlabel("-Log10 Adjusted P-value")
plt.ylabel("GO Terms")
plt.show()
