# The Stemness differentiation gene list was curated from https://doi.org/10.1038/ng.127
# The mucinous differentiation gene list was curated from https://doi.org/10.1016/j.cell.2021.11.017
# The mucinous differentiation gene list was curated from https://doi.org/10.1158/2159-8290.CD-24-0421


### Geneset analysis ###
## Stemness, Classical vs Basal, Mucinous differentiation ##

## Data preparation ##
f1 = "DEG_TCC_cal_4G.txt"
f2 = "sample_data.txt"
f3 = "Genelist.txt"

master_data = read.table(f1, header=TRUE, sep="\t", quote="", row.names = 1)
ccp = round(master_data[ ,1:12],digits = 0)
sample_names = colnames(ccp)
Gene_symbols = list(sub("^[^_]*_", "", rownames(ccp)))
names(Gene_symbols) = "Gene_symbols"

sample_names = c("AsPC1_Cont_1","AsPC1_Cont_2","AsPC1_Cont_3",
                 "AsPC1_Dox.1","AsPC1_Dox.2","AsPC1_Dox.3",
                 "AsPC1_MR_DOX_neg_1","AsPC1_MR_DOX_neg_2","AsPC1_MR_DOX_neg_3",
                 "AsPC1_MR_DOX_pos_1","AsPC1_MR_DOX_pos_2","AsPC1_MR_DOX_pos_3")
names_factor = c(1,1,1,2,2,2,3,3,3,4,4,4)

Sample_info = read.table(file = f2, header = TRUE, sep = "\t")
rownames(Sample_info) = Sample_info$Sample_names
ann_colors = list(
  Treatment = c(Cont = "#1F78B4", MUC1_KD = "#FB8072"),
  Status = c("Naïve" = "#66A61E", Resistant = "#7570B3")
)

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
  gene_list = Geneset_data[[i]]
  gene_list = gene_list[gene_list != ""]
  message("[", i, "] n_genes = ", length(gene_list))
  
  pre_data = ex_data[ex_data$Gene_symbols %in% gene_list, ]
  HM_data = pre_data[, 2:13, drop = FALSE]
  rownames(HM_data) = pre_data$Gene_symbols
  
  # Z-score scaling & clipping
  exccp_z = genescale(HM_data, axis = 1, method = "Z")
  exccp_z[exccp_z < -2] = -2
  exccp_z[exccp_z >  2] =  2
  
  # ==== Automatic sizing ====
  n_rows <- nrow(exccp_z)
  n_cols <- ncol(exccp_z)
  
  per_row_in  <- 0.11
  height_in   <- max(6, min(2.0 + n_rows * per_row_in, 30))
  
  top_keep_in    <- 1.1   # minimum top space
  bottom_keep_in <- 1.0   # minimum bottom space
  usable_in      <- max(2.0, height_in - top_keep_in - bottom_keep_in)  # height available for the heatmap body
  usable_mm      <- usable_in * 25.4
  
  # Adjust per-row maximum height (mm) to fit the available area
  row_mm_cap   <- 3.0                 # upper bound for row pitch
  row_mm_fit   <- usable_mm / n_rows  # row pitch derived from available area
  row_mm       <- min(row_mm_cap, row_mm_fit)
  
  # Use smaller fonts for row/column names as row count increases
  row_fs <- if (n_rows <= 60) 8 else if (n_rows <= 120) 7 else 6
  col_fs <- if (n_cols <= 12) 8 else 7
  
  # Maximum width of row names
  rn_max_w <- max_text_width(rownames(exccp_z), gp = gpar(fontsize = row_fs))
  
  # Fix top annotation height (prevent squashing)
  Annotation = HeatmapAnnotation(
    Treatment = Sample_info$Treatment,
    Status    = Sample_info$Status,
    show_annotation_name = TRUE,
    annotation_name_side = "right",
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
    annotation_legend_param = list(
      title_position= "topleft",
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 10)
    ),
    simple_anno_size = unit(3, "mm"),
    annotation_height = unit.c(unit(4, "mm"), unit(4, "mm")),
    col = ann_colors
  )
  
  col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  lgd = Legend(col_fun = col_fun, 
               title = "Z-scored\nexpression\nlevel", title_position = "topleft",
               title_gp = gpar(fontsize = 10, fontface = "bold"),
               labels = c("Low(-2)","High(2)"), labels_gp = gpar(fontsize = 10),
               at = c(-2, 2), grid_width = unit(3.5, "mm"), legend_height = unit(1.5, "cm")
  )
  
  # Heatmap body height = [#rows × row pitch]
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
    height = unit(n_rows * row_mm, "mm")
  )
  
  # TIFF output: device height is height_in (unchanged)
  tiff(file = paste0(i, ".tiff"), width = 10, height = height_in, units = "in", res = 600)
  set.seed(2024)
  draw(ht_with_names, padding = unit(c(4, 6, 6, 6), "mm"))  # Reserve a bit of outer padding
  draw(lgd, x = unit(0.982, "npc"), y = unit(0.87, "npc"), just = c("right", "top"))
  dev.off()
}



### Ex Basal and Classical gene set  ###
Geneset_data = read.table(f3, header=TRUE, sep="\t", quote="")
Basal_list = Geneset_data$scBasal[Geneset_data$scBasal != ""]
Classical_list = Geneset_data$scClassical[Geneset_data$scClassical != ""]

ex_data = cbind(Gene_symbols, ccp)

Basal_data = ex_data[ex_data$Gene_symbols %in% Basal_list, ]
Classical_data = ex_data[ex_data$Gene_symbols %in% Classical_list, ]

## Data preparation ##
Basal_HM_data = Basal_data[,2:13]
rownames(Basal_HM_data) = Basal_data$Gene_symbols
Classical_HM_data = Classical_data[,2:13]
rownames(Classical_HM_data) = Classical_data$Gene_symbols

# Align the common column (sample) order
common_cols <- intersect(colnames(Basal_HM_data), colnames(Classical_HM_data))
Basal_HM_data     <- as.matrix(Basal_HM_data[, common_cols, drop = FALSE])
Classical_HM_data <- as.matrix(Classical_HM_data[, common_cols, drop = FALSE])

if (exists("names_factor")) {
  if (!is.null(names(names_factor))) {
    names_factor <- names_factor[common_cols]
  }
}

library(genefilter)
z_scale_clip <- function(mat) {
  m <- genescale(mat, axis = 1, method = "Z")
  m[m < -2] <- -2
  m[m >  2] <-  2
  m
}
basal_z     <- z_scale_clip(Basal_HM_data)
classical_z <- z_scale_clip(Classical_HM_data)

library(circlize)
ann_colors = list(
  Treatment = c(Cont = "#1F78B4", MUC1_KD = "#FB8072"),
  Status = c("Naïve" = "#66A61E", Resistant = "#7570B3")
)

library(grid)
library(ComplexHeatmap)

Annotation = HeatmapAnnotation(
  Treatment = Sample_info$Treatment[match(common_cols, rownames(Sample_info))],
  Status    = Sample_info$Status[match(common_cols, rownames(Sample_info))],
  show_annotation_name = TRUE,
  annotation_name_side = "right",
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
  annotation_legend_param = list(
    title_position = "topleft",
    title_gp   = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 10)
  ),
  simple_anno_size = unit(3, "mm"),
  col = ann_colors
)

set.seed(2024)
col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

### Top panel: Basal
ht_basal <- Heatmap(
  basal_z,
  name = "Basal",
  row_title = "Basal gene set",
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  column_title = "Basal vs. Classical gene sets",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  column_order = colnames(basal_z),
  column_split = if (exists("names_factor")) names_factor else NULL,
  cluster_column_slices = FALSE,
  cluster_rows = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  column_names_side = "bottom",
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  col = col_fun,
  show_heatmap_legend = FALSE,
  row_dend_reorder = FALSE
)

### Bottom panel: Classical
ht_classical <- Heatmap(
  classical_z,
  name = "Classical",
  row_title = "Classical gene set",
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  column_title = "Basal vs. Classical gene sets",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  column_order = colnames(classical_z),
  column_split = if (exists("names_factor")) names_factor else NULL,
  cluster_column_slices = FALSE,
  cluster_rows = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  column_names_side = "bottom",
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  col = col_fun,
  show_heatmap_legend = FALSE,
  top_annotation = Annotation,
  row_dend_reorder = FALSE
)

# Output
tiff(file = "Basal_Classical_2row.tiff", width = 10, height = 8, units = "in", res = 600)
set.seed(2024)
ht_list <- ht_classical %v% ht_basal
draw(ht_list)
draw(lgd, x = unit(0.982, "npc"), y = unit(0.87, "npc"), just = c("right", "top"))
dev.off()

