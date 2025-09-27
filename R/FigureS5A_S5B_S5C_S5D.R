# The IFN-alpha/gamma gene list was curated from MgisgDB HALLMARK


### Extract IFN-alpha or IFN-gamma gene data ###
f1 = "DEG_TCC_cal_4G.txt"
f2 = "sample_data.txt"
f3 = "IFN_GENE_list.txt"

o1 = "Ex_IFN_alpha_in_M3.txt"
o2 = "Ex_IFN_gamma_in_M3.txt"
o3 = "Ex_IFN_alpha_in_M3_HM.tiff"
o4 = "Ex_IFN_gamma_in_M3_HM.tiff"
o5 = "Ex_IFN_alpha_in_M3_VP.tiff"
o6 = "Ex_IFN_gamma_in_M3_VP.tiff"

### Extract M3 genes ###
master_data = read.table(f1, header=TRUE, sep="\t", quote="", row.names = 1)
ex_tmp = subset(x = master_data, subset = estimatedDEG == 1)
ccp = round(ex_tmp[, 1:12], digits = 0)
sample_names = colnames(ccp)

### K-means clustering of genes ###
ex_ccp = ccp
library(genefilter)
exccp_z = genescale(ex_ccp, axis = 1, method = "Z")
set.seed(2024)
row_hc = kmeans(exccp_z, centers = 6)

library(tidyverse)
kmean_cluster_info = as.data.frame(row_hc$cluster)
kmean_gene = rownames(kmean_cluster_info)
row_hc_modi = row_hc$cluster %>%
  str_replace_all(c("1"="M1","2"="M2","3"="M3","4"="M4","5"="M5","6"="M6"))
names(row_hc_modi) = kmean_gene
alin_row_hc = factor(row_hc_modi, levels=c("M1","M2","M3","M4","M5","M6"))

# Extract gene list for module M3
ex_module = 3
Module_gene_list = names(row_hc$cluster[row_hc$cluster == ex_module])
length(Module_gene_list)
Module_data = ccp[Module_gene_list, ]

# Extract gene symbols from row names
Gene_symbols = list(sub("^[^_]*_", "", rownames(Module_data)))
names(Gene_symbols) = "Gene_symbols"


### Extract IFN-alpha or IFN-gamma gene sets ###
IFN_data = read.table(f3, header=TRUE, sep="\t", quote="")
IFN_ALPHA_list = IFN_data$GENE_ALPHA[IFN_data$GENE_ALPHA != ""]
IFN_GAMMA_list = IFN_data$GENE_GAMMA[IFN_data$GENE_GAMMA != ""]

# Combine symbols and expression for filtering
ex_data = cbind(Gene_symbols, Module_data)

alpha_data = ex_data[ex_data$Gene_symbols %in% IFN_ALPHA_list, ]
gamma_data = ex_data[ex_data$Gene_symbols %in% IFN_GAMMA_list, ]
write.table(alpha_data, file = o1, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)
write.table(gamma_data, file = o2, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)

## Significance testing via a meta-score using PCA (PC1) ##
set.seed(2025)
options(digits=3)
PCA_Pvalue = data.frame(row.names = c("2 <-> 1","3 <-> 1","4 <-> 1","3 <-> 2","4 <-> 2","4 <-> 3"))

# NOTE: define HM_alpha_data / HM_gamma_data before this block if running top-to-bottom
test_data_list = list(alpha = HM_alpha_data, gamma = HM_gamma_data)
library(multcomp)

for (i in seq_along(test_data_list)) {
  print(names(test_data_list[i]))
  ex_test_data = test_data_list[[i]]
  ex_test_data_z = genescale(ex_test_data, axis = 1, method = "Z")
  exdata_PCA = prcomp(ex_test_data_z, center = FALSE, scale. = FALSE)
  exdata_PC1 = data.frame(
    PC1 = exdata_PCA$rotation[, 1],
    cluster = as.factor(c(1,1,1,2,2,2,3,3,3,4,4,4))
  )
  fit = glm(PC1 ~ cluster, family = gaussian, data = exdata_PC1)
  multicomp = glht(fit, linfct = mcp(cluster = "Tukey"))
  test = summary(multicomp)
  p_value = ifelse(test$test$pvalues < 1e-16, "<1.00e-16", sprintf('%.3e', test$test$pvalues))
  coef = sprintf('%.2f', test$test$coefficients)
  a = as.data.frame(cbind(coef, p_value))
  colnames(a) = c(paste0(names(test_data_list[i]), "_coef"),
                  paste0(names(test_data_list[i]), "_Pvalue"))
  PCA_Pvalue = cbind(PCA_Pvalue, a)
}

comp_module = rownames(PCA_Pvalue)
write.table(cbind(comp_module, PCA_Pvalue),
            file = "IFN_PCA_Pvalue.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)



## Draw heatmap: IFN-alpha ##
## Data preparation##
HM_alpha_data = alpha_data[, 2:13]
rownames(HM_alpha_data) = alpha_data$Gene_symbols
library(genefilter)
exccp_z = genescale(HM_alpha_data, axis = 1, method = "Z")
exccp_z[exccp_z < -2] = -2
exccp_z[exccp_z >  2] =  2

Sample_info = read.table(file = f2, header = TRUE, sep = "\t")
rownames(Sample_info) = Sample_info$Sample_names

library(circlize)
ann_colors = list(
  Treatment = c(Cont = "#1F78B4", MUC1_KD = "#FB8072"),
  Status = c(Naïve = "#66A61E", Resistant = "#7570B3")
)

library(ggplot2)
library(ComplexHeatmap)

Annotation = HeatmapAnnotation(
  Treatment = Sample_info$Treatment,
  Status = Sample_info$Status,
  show_annotation_name = TRUE,
  annotation_name_side = "right",
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
  annotation_legend_param = list(
    title_position = "topleft",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 10)
  ),
  simple_anno_size = unit(3, "mm"),
  col = ann_colors
)

set.seed(2024)
sample_names = c("AsPC1_Cont_1","AsPC1_Cont_2","AsPC1_Cont_3",
                 "AsPC1_Dox.1","AsPC1_Dox.2","AsPC1_Dox.3",
                 "AsPC1_MR_DOX_neg_1","AsPC1_MR_DOX_neg_2","AsPC1_MR_DOX_neg_3",
                 "AsPC1_MR_DOX_pos_1","AsPC1_MR_DOX_pos_2","AsPC1_MR_DOX_pos_3")
names_factor = c(1,1,1,2,2,2,3,3,3,4,4,4)

gg = ComplexHeatmap::Heatmap(
  matrix = exccp_z,
  column_title = "Comparing IFN alpha genes in M3",
  column_order = colnames(exccp_z),
  cluster_column_slices = FALSE,
  column_split = names_factor,
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  column_names_side = "bottom",
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  heatmap_legend_param = list(
    at = c(-2, 2),
    labels = c("Low(-2)", "High(2)"),
    title = "Z-Scored Expression",
    title_position = "leftcenter-rot",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 10),
    legend_gp = list(x = unit(1, "npc"), y = unit(0.5, "npc"), just = c("center"))
  ),
  show_heatmap_legend = FALSE,
  top_annotation = Annotation,
  row_dend_reorder = FALSE
)

col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
lgd = Legend(
  col_fun = col_fun,
  title = "Z-Scored\nExpression\nLevel", title_position = "topleft",
  title_gp = gpar(fontsize = 10, fontface = "bold"),
  labels = c("Low(-2)", "High(2)"), labels_gp = gpar(fontsize = 10),
  at = c(-2, 2), grid_width = unit(3.5, "mm"), legend_height = unit(1.5, "cm")
)

tiff(file = o3, width = 10, height = 6, units = "in", res = 600)
set.seed(2024)
draw(gg)
draw(lgd, x = unit(0.982, "npc"), y = unit(0.87, "npc"), just = c("right", "top"))
dev.off()


## Violin plot with P-values (IFN-alpha) ##
set.seed(2025)
options(digits=3)
HM_alpha_data = alpha_data[, 2:13]
exccp_z = genescale(HM_alpha_data, axis = 1, method = "Z")

exdata_PCA = prcomp(exccp_z, center = FALSE, scale. = FALSE)
exdata_PC1 = data.frame(
  PC1 = exdata_PCA$rotation[, 1],
  cluster = as.factor(c(1,1,1,2,2,2,3,3,3,4,4,4))
)

aggregate(PC1 ~ cluster, data = exdata_PC1, mean)

fit = glm(PC1 ~ cluster, family = gaussian, data = exdata_PC1)
multicomp = glht(fit, linfct = mcp(cluster = "Tukey"))

M3_mean = as.data.frame(colMeans(exccp_z))
M3_mean$group = c("Cont","Cont","Cont","Dox","Dox","Dox",
                  "MR_DOX_neg","MR_DOX_neg","MR_DOX_neg",
                  "MR_DOX_pos","MR_DOX_pos","MR_DOX_pos")
y_max = round(max(M3_mean$`colMeans(exccp_z)`), 1)
y_min = round(min(M3_mean$`colMeans(exccp_z)`), 1)
p_high = abs(y_max - y_min) / 10
multcomp_p_df = data.frame(
  p_value = as.numeric(ifelse(summary(multicomp)$test$pvalues < 1e-16,
                              sprintf('%.3e', 1.00e-16),
                              sprintf('%.3e', summary(multicomp)$test$pvalues))),
  left_tip = c(1,1,1,2,2,3),
  right_tip = c(2,3,4,3,4,4),
  y_bar = c(y_max+p_high*1, y_max+p_high*2, y_max+p_high*3,
            y_max+p_high*4, y_max+p_high*5, y_max+p_high*6)
)

library(ggpmisc)
library(ggpp)
gg = ggplot(M3_mean, aes(x = M3_mean$group,
                         y = M3_mean$`colMeans(exccp_z)`,
                         fill = M3_mean$group)) +
  geom_violin(trim = FALSE) +
  scale_fill_manual(values = c("#f6aa00","#ff8082","#ff4b00","#804000"),
                    labels = c("1","2","3","4"), name = "") +
  geom_boxplot(width = .1, fill = "white") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x  = element_text(size = 12, face = "bold", colour = "black"),
    axis.title.y = element_text(size = 16),
    axis.text.y  = element_text(size = 14, face = "bold", colour = "black")
  ) +
  labs(title = "IFN alpha mean Z score on M3 in AsPC1 cell line",
       x = "Group", y = "Mean Z score") +
  scale_x_discrete(labels = c(
    "Cont" = "MRTX1133\nNaïve\n+ DMSO",
    "Dox"  = "MRTX1133\nNaïve\n+ DOX",
    "MR_DOX_neg" = "MRTX1133\nResistant\n+ DMSO",
    "MR_DOX_pos" = "MRTX1133\nResistant\n+ DOX"
  )) +
  geom_text_pairwise(
    data = multcomp_p_df,
    aes(xmin = left_tip, xmax = right_tip, y = y_bar,
        label = ifelse(p_value < 1.001e-16,
                       sprintf("P~`<`~%s", p_value),
                       sprintf("P~`=`~%s", p_value))),
    parse = TRUE, na.rm = TRUE, size = 4
  ) +
  expand_limits(y = max(M3_mean$`colMeans(exccp_z)`) + p_high)

ggsave(filename = o5, plot = gg, dpi = 600, width = 8, height = 6)
dev.off()



## IFN-gamma ##
HM_gamma_data = gamma_data[, 2:13]
rownames(HM_gamma_data) = gamma_data$Gene_symbols
library(genefilter)
exccp_z = genescale(HM_gamma_data, axis = 1, method = "Z")
exccp_z[exccp_z < -2] = -2
exccp_z[exccp_z >  2] =  2

gg = ComplexHeatmap::Heatmap(
  matrix = exccp_z,
  column_title = "Comparing IFN gamma genes in M3",
  column_order = colnames(exccp_z),
  cluster_column_slices = FALSE,
  column_split = names_factor,
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  column_names_side = "bottom",
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  heatmap_legend_param = list(
    at = c(-2, 2),
    labels = c("Low(-2)", "High(2)"),
    title = "Z-Scored Expression",
    title_position = "leftcenter-rot",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 10),
    legend_gp = list(x = unit(1, "npc"), y = unit(0.5, "npc"), just = c("center"))
  ),
  show_heatmap_legend = FALSE,
  top_annotation = Annotation,
  row_dend_reorder = FALSE
)

col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
lgd = Legend(
  col_fun = col_fun,
  title = "Z-Scored\nExpression\nLevel", title_position = "topleft",
  title_gp = gpar(fontsize = 10, fontface = "bold"),
  labels = c("Low(-2)", "High(2)"), labels_gp = gpar(fontsize = 10),
  at = c(-2, 2), grid_width = unit(3.5, "mm"), legend_height = unit(1.5, "cm")
)

tiff(file = o4, width = 10, height = 6, units = "in", res = 600)
set.seed(2024)
draw(gg)
draw(lgd, x = unit(0.982, "npc"), y = unit(0.87, "npc"), just = c("right", "top"))
dev.off()

## Violin plot with P-values (IFN-gamma) ##
set.seed(2025)
options(digits=3)
HM_gamma_data = gamma_data[, 2:13]
rownames(HM_gamma_data) = gamma_data$Gene_symbols
exccp_z = genescale(HM_gamma_data, axis = 1, method = "Z")

exdata_PCA = prcomp(exccp_z, center = FALSE, scale. = FALSE)
exdata_PC1 = data.frame(
  PC1 = exdata_PCA$rotation[, 1],
  cluster = as.factor(c(1,1,1,2,2,2,3,3,3,4,4,4))
)

aggregate(PC1 ~ cluster, data = exdata_PC1, mean)

fit = glm(PC1 ~ cluster, family = gaussian, data = exdata_PC1)
multicomp = glht(fit, linfct = mcp(cluster = "Tukey"))

M3_mean = as.data.frame(colMeans(exccp_z))
M3_mean$group = c("Cont","Cont","Cont","Dox","Dox","Dox",
                  "MR_DOX_neg","MR_DOX_neg","MR_DOX_neg",
                  "MR_DOX_pos","MR_DOX_pos","MR_DOX_pos")
y_max = round(max(M3_mean$`colMeans(exccp_z)`), 1)
y_min = round(min(M3_mean$`colMeans(exccp_z)`), 1)
p_high = abs(y_max - y_min) / 10
multcomp_p_df = data.frame(
  p_value = as.numeric(ifelse(summary(multicomp)$test$pvalues < 1e-16,
                              sprintf('%.3e', 1.00e-16),
                              sprintf('%.3e', summary(multicomp)$test$pvalues))),
  left_tip = c(1,1,1,2,2,3),
  right_tip = c(2,3,4,3,4,4),
  y_bar = c(y_max+p_high*1, y_max+p_high*2, y_max+p_high*3,
            y_max+p_high*4, y_max+p_high*5, y_max+p_high*6)
)

library(ggpmisc)
library(ggpp)
gg = ggplot(M3_mean, aes(x = M3_mean$group,
                         y = M3_mean$`colMeans(exccp_z)`,
                         fill = M3_mean$group)) +
  geom_violin(trim = FALSE) +
  scale_fill_manual(values = c("#f6aa00","#ff8082","#ff4b00","#804000"),
                    labels = c("1","2","3","4"), name = "") +
  geom_boxplot(width = .1, fill = "white") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x  = element_text(size = 12, face = "bold", colour = "black"),
    axis.title.y = element_text(size = 16),
    axis.text.y  = element_text(size = 14, face = "bold", colour = "black")
  ) +
  labs(title = "IFN gamma mean Z score on M3 in AsPC1 cell line",
       x = "Group", y = "Mean Z score") +
  scale_x_discrete(labels = c(
    "Cont" = "MRTX1133\nNaïve\n+ DMSO",
    "Dox"  = "MRTX1133\nNaïve\n+ DOX",
    "MR_DOX_neg" = "MRTX1133\nResistant\n+ DMSO",
    "MR_DOX_pos" = "MRTX1133\nResistant\n+ DOX"
  )) +
  geom_text_pairwise(
    data = multcomp_p_df,
    aes(xmin = left_tip, xmax = right_tip, y = y_bar,
        label = ifelse(p_value < 1.001e-16,
                       sprintf("P~`<`~%s", p_value),
                       sprintf("P~`=`~%s", p_value))),
    parse = TRUE, na.rm = TRUE, size = 4
  ) +
  expand_limits(y = max(M3_mean$`colMeans(exccp_z)`) + p_high)

ggsave(filename = o6, plot = gg, dpi = 600, width = 8, height = 6)
dev.off()

