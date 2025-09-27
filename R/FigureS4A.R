# Raw count data were generated with the STAR–RSEM pipeline from FASTQ files (GSE201412).
# Differentially expressed genes were identified using an edgeR-based procedure
# as implemented in multiple_comparison_test.R.
# The mucinous differentiation gene list was curated from https://doi.org/10.1158/2159-8290.CD-24-0421

### Analysis of mucinous differentiation ###
## Data preparation ##
f1 = "DEG_TCC_cal_vitro.txt"
f2 = "sample_data_vitro.txt"
f3 = "Genelist.txt"

master = read.table(f1, header=TRUE, sep="\t", quote="", row.names = 1)
log2_cpm = log2(master + 1)
ex_log2_cpm = log2_cpm[rowMeans(log2_cpm) > 0.32, ]

ccp = round(ex_log2_cpm[, 1:6], digits = 0)
sample_names = colnames(ccp)
Gene_symbols = list(sub("^[^_]*_", "", rownames(ccp)))
names(Gene_symbols) = "Gene_symbols"

names_factor = c(1,1,1,2,2,2)

Sample_info = read.table(file = f2, header = TRUE, sep = "\t")
rownames(Sample_info) = Sample_info$Sample_names

### Check stemness and mucinous differentiation gene sets ###
Geneset_data = read.table(f3, header=TRUE, sep="\t", quote="")
gene_set_list = colnames(Geneset_data)[c(1, 2, 5)]
ex_data = cbind(Gene_symbols, ccp)

library(genefilter)
library(circlize)
library(ggplot2)
library(ComplexHeatmap)
library(grid)

for (i in gene_set_list) {
  gene_list = get(paste0(i, "_genes"))
  gene_list = gene_list[gene_list != ""]
  message("[", i, "] n_genes = ", length(gene_list))
  
  pre_data = ex_data[ex_data$Gene_symbols %in% gene_list, ]
  HM_data = pre_data[, 2:6, drop = FALSE]
  rownames(HM_data) = pre_data$Gene_symbols
  
  # Z-score scaling and clipping
  exccp_z = genescale(HM_data, axis = 1, method = "Z")
  exccp_z[exccp_z < -2] = -2
  exccp_z[exccp_z >  2] =  2
  
  # ==== Automatic sizing (with row names shown) ====
  n_rows <- nrow(exccp_z)
  n_cols <- ncol(exccp_z)
  
  # Device height (with row names) scales with the number of rows; cap at 30 inches
  per_row_in  <- 0.11
  height_in   <- max(6, min(2.0 + n_rows * per_row_in, 30))
  
  # Reserve vertical space for title/legend/annotations/column names (top/bottom)
  top_keep_in    <- 1.1   # minimum space at the top
  bottom_keep_in <- 1.0   # minimum space at the bottom
  usable_in      <- max(2.0, height_in - top_keep_in - bottom_keep_in)  # height available for the heatmap body
  usable_mm      <- usable_in * 25.4
  
  # Per-row height (mm), fitted to the available area
  row_mm_cap   <- 3.0                 # upper bound for row pitch
  row_mm_fit   <- usable_mm / n_rows  # data-driven row pitch
  row_mm       <- min(row_mm_cap, row_mm_fit)
  
  # Font sizes for row/column names adjusted by table size
  row_fs <- if (n_rows <= 60) 8 else if (n_rows <= 120) 7 else 6
  col_fs <- if (n_cols <= 12) 8 else 7
  
  # Maximum width required for row names
  rn_max_w <- max_text_width(rownames(exccp_z), gp = gpar(fontsize = row_fs))
  
  # Fix top annotation height (prevent squashing)
  Annotation = HeatmapAnnotation(
    Treatment = Sample_info$Treatment,
    Status    = Sample_info$Status,
    show_annotation_name = TRUE,
    annotation_name_side = "right",
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
    annotation_legend_param = list(
      title_position = "topleft",
      title_gp  = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 10)
    ),
    simple_anno_size = unit(3, "mm"),
    annotation_height = unit.c(unit(4, "mm"), unit(4, "mm")),
    col = ann_colors
  )
  
  # Heatmap body height = [n_rows × row pitch]
  ht_with_names <- Heatmap(
    matrix = exccp_z,
    name = i,
    column_title = paste0("Comparing ", i),
    column_order = colnames(exccp_z),
    cluster_column_slices = FALSE,
    column_split = names_factor,
    cluster_rows = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "ward.D2",
    column_names_side = "bottom",
    column_names_gp = gpar(fontsize = col_fs),
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = row_fs),
    row_names_max_width = rn_max_w,
    col = col_fun,
    show_heatmap_legend = FALSE,
    top_annotation = Annotation,
    row_dend_reorder = FALSE,
    height = unit(n_rows * row_mm, "mm")  # fits within the available area
  )
  
  # TIFF output: device height is height_in (unchanged)
  tiff(file = paste0(i, ".tiff"), width = 10, height = height_in, units = "in", res = 600)
  set.seed(2024)
  draw(ht_with_names, padding = unit(c(4, 6, 6, 6), "mm"))
  draw(lgd, x = unit(0.982, "npc"), y = unit(0.87, "npc"), just = c("right", "top"))
  dev.off()
}
