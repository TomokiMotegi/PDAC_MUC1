### ConsensusClusterPlus on samples ###
f1 = "DEG_TCC_cal_4G.txt"
f2 = "sample_data.txt"
f3 = "4_Group_Clustering"

master_data = read.table(f1, header = TRUE, sep = "\t", quote = "", row.names = 1)
ex_tmp = subset(x = master_data, subset = estimatedDEG == 1)

# Use the first 12 sample columns as the expression matrix
ccp = round(ex_tmp[, 1:12], digits = 0)
sample_names = colnames(ccp)

library(ConsensusClusterPlus)
rcc = ConsensusClusterPlus(
  as.matrix(ccp),
  maxK = 5,
  reps = 100,
  pItem = 0.8,
  pFeature = 1,
  seed = 2024,
  distance = "euclidean",
  clusterAlg = "km",
  innerLinkage = "ward.D2",  # used for hierarchical methods; harmless if unused
  finalLinkage = "ward.D2",  # used for hierarchical methods; harmless if unused
  plot = "png",
  title = f3
)
# Select the desired number of clusters to inspect
nc = 4

### K-means clustering of DEGs ###
ex_ccp = ccp
library(genefilter)
exccp_z = genescale(ex_ccp, axis = 1, method = "Z")
set.seed(2024)
row_hc = kmeans(exccp_z, centers = 6)

library(tidyverse)
kmean_cluster_info = as.data.frame(row_hc$cluster)
kmean_gene = rownames(kmean_cluster_info)

# Rename clusters to M1..M6 for readability
row_hc_modi = row_hc$cluster %>%
  str_replace_all(c("1" = "M1", "2" = "M2", "3" = "M3",
                    "4" = "M4", "5" = "M5", "6" = "M6"))
names(row_hc_modi) = kmean_gene
alin_row_hc = factor(row_hc_modi, levels = c("M1", "M2", "M3", "M4", "M5", "M6"))

## Sample-side metadata ##
Sample_info = read.table(file = f2, header = TRUE, sep = "\t")

# Append consensus class (sample clusters) for the chosen K
cluster = as.data.frame(rcc[[nc]]$consensusClass)
colnames(cluster) = "cluster"
Sample_info = cbind(Sample_info, cluster)
rownames(Sample_info) = Sample_info$Sample_names


### G3 vs G4 comparison ###
add_G3_G4_data <- ccp %>%
  mutate(
    Mean_MR_DOX_neg = rowMeans(across(starts_with("AsPC1_MR_DOX_neg"))),
    Mean_MR_DOX_pos = rowMeans(across(starts_with("AsPC1_MR_DOX_pos"))),
    Fold_Change = ifelse(Mean_MR_DOX_neg == 0, NA, Mean_MR_DOX_pos / (Mean_MR_DOX_neg + 1e-6)),
    m.value = ifelse(is.na(Fold_Change), NA, log2(Fold_Change))
  )

# Extract gene IDs after the first underscore in row names
Ensembl_ID = sub("^[^_]*_", "", rownames(add_G3_G4_data))
analysis_data = cbind(Ensembl_ID, add_G3_G4_data$m.value)
rownames(analysis_data) = rownames(add_G3_G4_data)

## Create a ranked vector for GSEA ##
library(tidyverse)

mapped_df = na.omit(as.data.frame(analysis_data))
colnames(mapped_df) = c("external_gene_name", "m.value")
mapped_df$m.value = as.numeric(mapped_df$m.value)
str(mapped_df$m.value)

# Check duplicated symbols
any(duplicated(mapped_df$external_gene_name))

dup_gene_symbols <- mapped_df %>%
  dplyr::filter(duplicated(external_gene_name)) %>%
  dplyr::pull(external_gene_name)

mapped_df %>%
  dplyr::filter(external_gene_name %in% dup_gene_symbols) %>%
  dplyr::arrange(external_gene_name)

filtered_mapped_df <- mapped_df %>%
  dplyr::arrange(dplyr::desc(abs(m.value))) %>%
  dplyr::distinct(external_gene_name, .keep_all = TRUE)

# Confirm no duplicates remain
any(duplicated(filtered_mapped_df$external_gene_name))

# Named vector ranked by log2 fold change
lfc_vector <- filtered_mapped_df$m.value
names(lfc_vector) <- filtered_mapped_df$external_gene_name
lfc_vector <- sort(lfc_vector, decreasing = TRUE)
head(lfc_vector)

any(!is.finite(lfc_vector))
lfc_vector <- lfc_vector[is.finite(lfc_vector)]


### MSigDB Hallmark GSEA ###
library(msigdbr)
MSigDB_category = msigdbr_collections()

MSigDB_cat = "Hallmark_G3vsG4"
Msig_data <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)

library(clusterProfiler)
set.seed(2023)

# Run GSEA against Hallmark gene sets using the ranked vector
gsea_results <- GSEA(
  geneList = lfc_vector,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  eps = 0,
  seed = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE = dplyr::select(Msig_data, gs_name, gene_symbol)
)

# Create an output folder tagged by the MSigDB category
if (!dir.exists(paste0(gsub(":", "_", MSigDB_cat)))) {
  dir.create(paste0(gsub(":", "_", MSigDB_cat)))
}

gsea_result_df = data.frame(gsea_results@result)


## Visualization and export ##
library(ggplot2)
library(enrichplot)

gsea_result_df %>%
  dplyr::slice_max(NES, n = 3)

## GSEA plots for selected Hallmark terms ##
name = "HALLMARK_KRAS_SIGNALING_UP"
ex_gene_set = subset(x = gsea_result_df, subset = ID == name)

most_positive_nes_plot <- enrichplot::gseaplot(
  gsea_results,
  geneSetID = name,
  title = name,
  color.line = "#0d76ff"
)

ggsave(
  filename = paste0("./", gsub(":", "_", MSigDB_cat), "/",
                    gsub(":", "_", MSigDB_cat), "_", name, ".tiff"),
  plot = most_positive_nes_plot, width = 6, height = 8, dpi = 600
)

name = "HALLMARK_KRAS_SIGNALING_DN"
ex_gene_set = subset(x = gsea_result_df, subset = ID == name)

most_positive_nes_plot <- enrichplot::gseaplot(
  gsea_results,
  geneSetID = name,
  title = name,
  color.line = "#0d76ff"
)

ggsave(
  filename = paste0("./", gsub(":", "_", MSigDB_cat), "/",
                    gsub(":", "_", MSigDB_cat), "_", name, ".tiff"),
  plot = most_positive_nes_plot, width = 6, height = 8, dpi = 600
)
