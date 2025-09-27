### Fold-change plots ###
f1 = "DEG_TCC_cal_4G.txt"
master_data = read.table(f1, header = TRUE, sep = "\t", quote = "", row.names = 1)
ex_tmp = subset(x = master_data, subset = estimatedDEG == 1)
ccp = round(ex_tmp[, 1:12], digits = 0)

# Compute per-condition means and log2 fold changes
mean_muc1kd_fc_data <- ccp %>%
  mutate(
    a.value = log2(rowMeans(x = ccp, na.rm = TRUE)),
    Mean_Cont = rowMeans(across(starts_with("AsPC1_Cont")), na.rm = TRUE),        # mean of control replicates
    Mean_DOX  = rowMeans(across(starts_with("AsPC1_Dox")),  na.rm = TRUE),        # mean of DOX-treated replicates
    Mean_MR_DOX_neg = rowMeans(across(starts_with("AsPC1_MR_DOX_neg")), na.rm = TRUE),  # mean of MR DOX− replicates
    Mean_MR_DOX_pos = rowMeans(across(starts_with("AsPC1_MR_DOX_pos")), na.rm = TRUE),  # mean of MR DOX+ replicates
    AsPC1_FC    = log2(ifelse(Mean_Cont == 0, NA, Mean_DOX / (Mean_Cont + 1e-6))),      # add small epsilon to avoid division by zero
    AsPC1_MR_FC = log2(ifelse(Mean_MR_DOX_neg == 0, NA, Mean_MR_DOX_pos / (Mean_MR_DOX_neg + 1e-6)))  # avoid division by zero
  )

# Extract gene identifiers from row names
gene_ids = sub("^[^_]*_", "", rownames(mean_muc1kd_fc_data))

# Assemble a tidy data frame for plotting
kras_fig_mean_expr_df <- data.frame(
  Genes       = sub("^[^_]*_", "", rownames(mean_muc1kd_fc_data)),
  Mean        = as.numeric(mean_muc1kd_fc_data$a.value),
  AsPC1_FC    = as.numeric(mean_muc1kd_fc_data$AsPC1_FC),
  AsPC1_MR_FC = as.numeric(mean_muc1kd_fc_data$AsPC1_MR_FC),
  row.names   = NULL
)

# Keep rows with finite FC values for both comparisons
kras_fig_mean_expr_df <- kras_fig_mean_expr_df %>%
  filter(is.finite(AsPC1_FC) & is.finite(AsPC1_MR_FC))

# Load MSigDB Hallmark gene sets
library(msigdbr)
Msig_data <- msigdbr(
  species  = "Homo sapiens",
  category = "H"
)

# KRAS Hallmark gene sets (UP/DN)
kras_up_geneset = subset(x = Msig_data, subset = Msig_data$gs_name == "HALLMARK_KRAS_SIGNALING_UP")
kras_up_symbol  = unique(kras_up_geneset$gene_symbol)
length(kras_up_symbol)

kras_dn_geneset = subset(x = Msig_data, subset = Msig_data$gs_name == "HALLMARK_KRAS_SIGNALING_DN")
kras_dn_symbol  = unique(kras_dn_geneset$gene_symbol)
length(kras_dn_symbol)

kras_symbol = unique(sort(c(kras_up_symbol, kras_dn_symbol)))

# Plotting libraries
library(ggplot2)
library(dplyr)
library(ggrepel)


## KRAS UP ##
# Subset to genes in HALLMARK_KRAS_SIGNALING_UP
ex_kras_fig = subset(x = kras_fig_mean_expr_df,
                     subset = kras_fig_mean_expr_df$Genes %in% kras_up_symbol)

# Color by concordant up/down vs other
ex_kras_fig <- ex_kras_fig %>%
  mutate(color_group = case_when(
    AsPC1_FC > 1  & AsPC1_MR_FC > 1  ~ "Up in Both",
    AsPC1_FC < -1 & AsPC1_MR_FC < -1 ~ "Down in Both",
    TRUE                              ~ "Other"
  ))

gg = ggplot(ex_kras_fig, aes(x = AsPC1_FC, y = AsPC1_MR_FC, color = color_group)) +
  geom_point(alpha = 0.6) +
  geom_text_repel(
    data = filter(ex_kras_fig, color_group %in% c("Up in Both", "Down in Both")),
    aes(label = paste0(Genes, "_", round(Mean, 2))),
    size = 3,
    box.padding = 0.4,
    max.overlaps = 30,
    show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(
    values = c("Up in Both" = "red", "Down in Both" = "blue", "Other" = "gray70"),
    name   = "",
    labels = c("Up in Both" = "|FC| ≥ 1 in both", "Down in Both" = "≤ −1 in both", "Other" = "Other")
  ) +
  labs(
    title = "Scatter plot of fold changes with gene labels and mean expression",
    x = "AsPC1 MUC1 KD FC (log2)",
    y = "AsPC1/MRTX1133 MUC1 KD FC (log2)"
  ) +
  theme_classic()

ggsave(filename = "kras_up_label.tiff", plot = gg, dpi = 600, width = 8, height = 8)
dev.off()


## KRAS DN ##
# Subset to genes in HALLMARK_KRAS_SIGNALING_DN
ex_kras_fig = subset(x = kras_fig_mean_expr_df,
                     subset = kras_fig_mean_expr_df$Genes %in% kras_dn_symbol)

# Color by concordant up/down vs other
ex_kras_fig <- ex_kras_fig %>%
  mutate(color_group = case_when(
    AsPC1_FC > 1  & AsPC1_MR_FC > 1  ~ "Up in Both",
    AsPC1_FC < -1 & AsPC1_MR_FC < -1 ~ "Down in Both",
    TRUE                              ~ "Other"
  ))

gg = ggplot(ex_kras_fig, aes(x = AsPC1_FC, y = AsPC1_MR_FC, color = color_group)) +
  geom_point(alpha = 0.6) +
  geom_text_repel(
    data = filter(ex_kras_fig, color_group %in% c("Up in Both", "Down in Both")),
    aes(label = paste0(Genes, "_", round(Mean, 2))),
    size = 3,
    box.padding = 0.4,
    max.overlaps = 30,
    show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(
    values = c("Up in Both" = "red", "Down in Both" = "blue", "Other" = "gray70"),
    name   = "",
    labels = c("Up in Both" = "|FC| ≥ 1 in both", "Down in Both" = "≤ −1 in both", "Other" = "Other")
  ) +
  labs(
    title = "HALLMARK_KRAS_SIGNALING_DN",
    x = "AsPC1 MUC1 KD FC (log2)",
    y = "AsPC1/MRTX1133 MUC1 KD FC (log2)"
  ) +
  theme_classic()

ggsave(filename = "kras_dn_label.tiff", plot = gg, dpi = 600, width = 8, height = 8)
