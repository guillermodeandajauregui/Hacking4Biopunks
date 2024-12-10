# -------------------------------------
# 1. Cargar librerías
# -------------------------------------
library(tidyverse)
library(vroom)
library(pheatmap)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)

# -------------------------------------
# 2. Preparación de datos
# -------------------------------------
# Cargar datos de expresión
expr_data <- vroom::vroom("/datos/rnaseq_demos/GSE67333_norm_counts_FPKM_GRCh38.p13_NCBI.tsv.gz")

# Cargar anotaciones de muestras
annotations <- read.csv("/datos/rnaseq_demos/gse67333_labels.txt", header = TRUE) %>%
  mutate(Condition = str_sub(Condition, end = -2))  # Limpiar etiquetas de condición

# Asegurar que las columnas coincidan con las anotaciones
expr_data <- expr_data %>%
  column_to_rownames(var = "GeneID") %>%
  select(all_of(annotations$SampleID))  # Ordenar según las anotaciones

# -------------------------------------
# 3. Visualización inicial: Heatmap
# -------------------------------------
# Subconjunto de genes más variables
top_genes <- expr_data %>%
  mutate(variance = apply(., 1, var)) %>%
  arrange(desc(variance)) %>%
  slice_head(n = 50) %>%
  select(-variance)

# Generar heatmap
pheatmap(
  mat = top_genes,
  show_rownames = FALSE,
  show_colnames = TRUE
)

# -------------------------------------
# 4. Análisis de PCA
# -------------------------------------
# Preparar datos para PCA
pca_data <- expr_data %>%
  t() %>%
  prcomp(center = TRUE, scale. = FALSE)

# Integrar PCA con anotaciones
pca_df <- as_tibble(pca_data$x) %>%
  mutate(SampleID = rownames(pca_data$x)) %>%
  left_join(annotations, by = "SampleID")

# Scatterplot del PCA
ggplot(pca_df) +
  aes(x = PC1, y = PC2, color = Condition, label = SampleID) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(nudge_y = 5, size = 3) +
  labs(
    title = "PCA de las muestras",
    x = "PC1",
    y = "PC2"
  ) +
  theme_minimal()

# -------------------------------------
# 5. Análisis de expresión diferencial con limma
# -------------------------------------
# Definir diseño del modelo
design <- model.matrix(~Condition, data = annotations)

# Ajuste del modelo lineal
fit <- lmFit(expr_data, design)
fit <- eBayes(fit)

# Resultados completos
limma_res <- topTable(fit, coef = "ConditionLOAD", number = Inf, sort.by = "none") %>%
  rownames_to_column("gene_id") %>%
  as_tibble()

# Filtro de genes significativos
limma_res_filtered <- limma_res %>%
  filter(B > 1, abs(logFC) > 1)

# Guardar resultados filtrados
vroom::vroom_write(limma_res_filtered, "/datos/rnaseq_demos/diffex_demo.txt")

# -------------------------------------
# 6. Enriquecimiento funcional con Gene Ontology
# -------------------------------------
# Leer anotaciones de genes
annot <- vroom::vroom("/datos/rnaseq_demos/Human.GRCh38.p13.annot.tsv.gz")

# Mapear IDs de genes a símbolos HGNC
diffex_with_hgnc <- limma_res_filtered %>%
  left_join(annot, by = c("gene_id" = "GeneID")) %>%
  select(gene_id, logFC, P.Value, adj.P.Val, Symbol) %>%
  filter(!is.na(Symbol))

# Crear lista de genes
gene_list <- diffex_with_hgnc$Symbol

# Enriquecimiento en GO
go_enrichment <- enrichGO(
  gene = gene_list,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH"
)

# Visualizar resultados
go_enrichment@result %>%
  as_tibble() %>%
  filter(qvalue < 0.001) %>%
  ggplot() +
  aes(x = -log10(qvalue), y = Description) +
  geom_bar(stat = "identity") +
  labs(title = "Enriquecimiento de GO", x = "-Log10 Q-value", y = "Procesos biológicos") +
  theme_minimal()
