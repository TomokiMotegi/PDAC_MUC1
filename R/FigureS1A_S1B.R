# Raw count data were generated with the STAR–RSEM pipeline from FASTQ files (GSE201412).
# Differentially expressed genes were identified using the same procedure as in multiple_comparison_test.R.

### Volcano plots highlighting MUC1 and ALPG ###
res_type = c("virto", "vivo")  # NOTE: if intended, "vitro" may be the correct spelling
library(ggplot2)
library(dplyr)
library(stringr)
library(tibble)
library(ggrepel)

for (i in res_type) {
  in_f = paste0("DEG_TCC_cal_",i,".txt")
  out_f = paste0("AsPC1_MRT24h_",i,"_Volcano.tiff")
  master = read.table(in_f, header=TRUE, row.names=1, check.names = TRUE, 
                      sep="\t", quote="") 
  tmp = master
  # Use log2 fold change and P-values for volcano plotting
  FC_z <- tmp$logFC
  vol_data <- tmp %>%
    rownames_to_column(var = "rownames") %>%
    select(P.Value, rownames) %>%
    mutate(FC_z = tmp$logFC)
  
  sig_data <- vol_data %>% filter(P.Value < 0.05 & (FC_z >= 1 | FC_z <= -1))
  # Target genes to annotate (MUC1, ALPG)
  target_data = subset(vol_data,rownames == "ENSG00000185499_MUC1" | 
                           rownames == "ENSG00000163286_ALPG")
  
  max_right_p <- max(-log10(vol_data %>% filter(FC_z >= 0) %>% pull(P.Value)))
  max_left_p <- max(-log10(vol_data %>% filter(FC_z < 0) %>% pull(P.Value)))
  max_fc <- max(vol_data$FC_z)
  min_fc <- min(vol_data$FC_z)
  
  # Volcano plot
  gg <- ggplot(vol_data, aes(x = FC_z, y = -log10(P.Value))) +
    geom_point(color = "lightgrey", size = 2) + 
    geom_point(data = sig_data %>% filter(FC_z < 1), 
               aes(x = FC_z, y = -log10(P.Value)), 
               color = "deepskyblue", size = 2) + 
    geom_point(data = sig_data %>% filter(FC_z > 1), 
               aes(x = FC_z, y = -log10(P.Value)), 
               color = "orange", size = 2) + 
    geom_point(data = cytokine_data, aes(x = FC_z, 
                                         y = -log10(P.Value)), 
               color = "red", size = 3) + 
    geom_text_repel(data = cytokine_data, 
                    aes(x = FC_z, y = -log10(P.Value), label = str_sub(rownames, 17, -1)),
                    size = 6, 
                    fontface = "bold", 
                    box.padding = 0.5, 
                    point.padding = 0.5, 
                    segment.color = "black",
                    segment.size = 0.5,
                    max.overlaps = Inf,
                    force = 10,
                    nudge_y = 0.2) +
    geom_vline(xintercept = c(1, -1), color = "grey", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), color = "grey", linetype = "dashed") +
    labs(
      title = "",
      x = "Log2 Fold Change",
      y = "-Log10(P Value)"
    ) +
    theme_classic() +
    theme(
      text = element_text(size = 14),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14, color = "black"),
      panel.grid = element_blank()
    )
  
  ggsave(filename = out_f, plot = gg, width = 6, height = 6, dpi = 600)
}

