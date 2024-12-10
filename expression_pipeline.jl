# -------------------------------------
# 1. Load Required Libraries
# -------------------------------------
using DataFrames
using CSV
using StatsBase
using Plots

# -------------------------------------
# 2. Load Data
# -------------------------------------
# Load expression data
expr_data = CSV.read("/datos/rnaseq_demos/GSE67333_norm_counts_FPKM_GRCh38.p13_NCBI.tsv.gz", DataFrame)
rename!(expr_data, :GeneID => :gene_id)  # Rename the first column

# Load sample annotations
annotations = CSV.read("/datos/rnaseq_demos/gse67333_labels.txt", DataFrame)
annotations.Condition = replace.(annotations.Condition, r"(\d+)" => "")

# Ensure columns match annotations
expr_data = select(expr_data, [:gene_id, annotations.SampleID...])

# -------------------------------------
# 3. Heatmap of Top Variable Genes
# -------------------------------------
# Calculate variance for each gene
expr_data.variance = [var(row[2:end]) for row in eachrow(expr_data)]

# Select top 50 variable genes
top_genes = expr_data |> sort(:variance, rev = true) |> first(50)

# Plot heatmap
heatmap(Matrix(select(top_genes, Not([:gene_id, :variance]))),
    title = "Top 50 Variable Genes Heatmap",
    xlabel = "Samples",
    ylabel = "Genes",
    color=:viridis)

# -------------------------------------
# 4. PCA Analysis
# -------------------------------------
# Transpose expression data for PCA
sample_matrix = Matrix(select(expr_data, Not([:gene_id, :variance])))'
gene_ids = expr_data.gene_id

# Perform PCA
using MultivariateStats
pca_model = fit(PCA, sample_matrix, maxoutdim = 2)
pca_scores = projection(pca_model)

# Merge PCA results with annotations
pca_df = DataFrame(PC1 = pca_scores[:, 1], PC2 = pca_scores[:, 2], SampleID = annotations.SampleID, Condition = annotations.Condition)

# Plot PCA scatterplot
scatter(pca_df.PC1, pca_df.PC2,
    group=pca_df.Condition,
    title="PCA of Samples",
    xlabel="PC1",
    ylabel="PC2",
    label="Condition",
    legend=:topright)

# -------------------------------------
# 5. Differential Expression Analysis
# -------------------------------------
# Define groups
group_ctrl = findall(x -> x == "CTRL", annotations.Condition)
group_load = findall(x -> x == "LOAD", annotations.Condition)

# Compute log2 fold change
logFC = [mean(log2.(expr_data[group_load, i] .+ 1)) - mean(log2.(expr_data[group_ctrl, i] .+ 1)) for i in 2:size(expr_data, 2)]

# Perform t-tests
using HypothesisTests
p_values = [pvalue(TTest(log2.(expr_data[group_load, i] .+ 1), log2.(expr_data[group_ctrl, i] .+ 1))) for i in 2:size(expr_data, 2)]

# Adjust p-values (Benjamini-Hochberg)
using FDR
adj_p_values = FDR.adjust(p_values)

# Create results DataFrame
diffex_results = DataFrame(
    gene_id = expr_data.gene_id,
    logFC = logFC,
    PValue = p_values,
    adjPValue = adj_p_values
)

# Filter significant genes
filtered = filter(row -> abs(row.logFC) > 1 && row.adjPValue < 0.05, diffex_results)

# Save filtered results
CSV.write("/datos/rnaseq_demos/diffex_demo_julia.txt", filtered)

# -------------------------------------
# 6. Functional Enrichment Analysis
# -------------------------------------
# Map genes to HGNC symbols (assumes annotation file)
annot = CSV.read("/datos/rnaseq_demos/Human.GRCh38.p13.annot.tsv.gz", DataFrame)

# Merge significant genes with annotations
filtered_annot = leftjoin(filtered, annot, on = :gene_id)

# Perform a basic overrepresentation test (hypergeometric test)
function hypergeometric_enrichment(gene_set, background_set, all_genes)
    N = length(all_genes)  # Total genes
    K = length(background_set)  # Background size
    n = length(gene_set)  # Observed size
    x = length(intersect(gene_set, background_set))  # Intersection size
    Hypergeometric(K, N - K, n)(x)
end

# Example GO background (replace with real data)
example_go_set = ["GeneA", "GeneB", "GeneC", "GeneD"]
all_genes = annot.Symbol

# Run enrichment
enrichment_results = hypergeometric_enrichment(filtered_annot.Symbol, example_go_set, all_genes)

# -------------------------------------
# 7. Visualize Enrichment Results
# -------------------------------------
# (Here we mock enrichment visualization, adjust as needed)
bar(["GO Term 1", "GO Term 2"], [0.01, 0.05],
    title="Enriched GO Terms",
    xlabel="GO Terms",
    ylabel="-log10(P-Value)")
