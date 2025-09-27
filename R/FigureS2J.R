# The list of differentially expressed genes (DEGs) was obtained by comparing
# AsPC1_MR_DOX_neg vs AsPC1_MR_DOX_pos using an edgeR-based pipeline
# implemented in multiple_comparison_test.R.

### GSEA for AsPC-1/MR/tet-MUC1shRNA cells treated with vehicle or DOX for 7 days ###
f1 = "AsPC1_MR_MUC1KD_DEG_TCC_cal.txt"
master_data = read.table(f1, header=TRUE, sep="\t", quote="")
head(master_data)

DEG_data = subset(master_data, a.value > 2)
Ensembl_ID = sub("^[^_]*_", "", DEG_data$rownames.tcc.count.)
analysis_data = cbind(Ensembl_ID, 
                      DEG_data[,9:14])

## Map to gene symbols and prepare a ranked list ##
library(tidyverse)
mapped_df = analysis_data
names(mapped_df)[1] = "external_gene_name"

# Check for duplicated gene symbols
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

# Re-check that no duplicates remain
any(duplicated(filtered_mapped_df$external_gene_name))


# Create a named vector ranked by the log2 fold change values
lfc_vector <- filtered_mapped_df$m.value
names(lfc_vector) <- filtered_mapped_df$external_gene_name
lfc_vector <- sort(lfc_vector, decreasing = TRUE)


### Load MSigDB collections and Hallmark gene sets ###
library(msigdbr)
MSigDB_cat = "AsPC1_MR_MUC1KD"
Msig_data <- msigdbr(
  species = "Homo sapiens",
  category = "H"             # Hallmark gene sets
)

library(clusterProfiler)
set.seed(2023)

# Run GSEA against the Hallmark gene sets
gsea_results <- GSEA(
  geneList = lfc_vector,
  minGSSize = 25, 
  maxGSSize = 500, 
  pvalueCutoff = 0.50, 
  eps = 0,  
  seed = TRUE, 
  pAdjustMethod = "BH",  
  TERM2GENE = dplyr::select(Msig_data, gs_name, gene_symbol)
)

# Create an output folder if it does not already exist
if(!dir.exists(paste0(gsub(":", "_", MSigDB_cat)))){
  dir.create(paste0(gsub(":", "_", MSigDB_cat)))
}

library(enrichplot)

NES_min = floor(min(gsea_results@result$NES))
NES_max = ceiling(max(gsea_results@result$NES))

# Dot plot of GSEA results (split by sign of NES)
p_go = dotplot(gsea_results, showCategory = 9, 
               split = ".sign", x = "NES", color = "p.adjust", 
               decreasing = TRUE, 
               font.size = 8, label_format = 30) + 
  scale_x_continuous(limit = c(NES_min, NES_max),  
                     breaks = seq(NES_min, NES_max, 1)) +
  scale_y_discrete(labels = function(labels) {
    str_remove(labels, "HALLMARK_") 
  }) + 
  geom_vline(xintercept = 0, linetype = "solid", 
             color = "#444444", size = 1) +
  scale_fill_gradientn(colors = c("red", "blue"),
                       name = "p.adjust") +
  theme( 
    legend.text = element_text(size = 10),
    legend.position = "bottom",
    legend.box = "vertical",
    axis.title.x = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    text = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

ggsave(filename = paste0("./", gsub(":", "_", MSigDB_cat), "/", 
                         gsub(":", "_", MSigDB_cat), "_qvalue_dotplot.tiff"), 
       width = 6, height = 6, dpi = 600)
