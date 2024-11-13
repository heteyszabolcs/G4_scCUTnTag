print("Load R packages")
if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(
  "data.table",
  "Seurat",
  "ComplexHeatmap",
  "circlize",
  "glue",
  "tidyverse",
  "EnhancedVolcano",
  "ggpubr",
  "ggrastr"
)

# helper function for annotation
source("../../utils/annotation.R")

# folders
result_folder = "../../results/Seurat/"

# logistic regression based differential analysis
# find G4 marker regions for cluster 0
g4 = readRDS("../../results/Seurat/unsorted_mousebrain/res0.1/outputs/Seurat_object.Rds")
markers_for_cl0 = FindMarkers(g4, test.use = "LR", latent.vars = "peak_region_fragments", ident.1 = 0, ident.2 = 1, group.by = "seurat_clusters") 
markers_for_cl0 = markers_for_cl0 %>% mutate(id = rownames(markers_for_cl0))
markers_for_cl0 = markers_for_cl0 %>% separate(id, sep = "-", into = c("chr", "start", "end"), remove = TRUE) %>% 
  mutate(start = as.numeric(start), end = as.numeric(end))
markers_for_cl0 = markers_for_cl0 %>% mutate(region = paste(chr, start, end, sep = "-"))
markers_for_cl0_annot = mm10_annotation(regions = markers_for_cl0, start_col = "start", end_col = "end", seqname_col = "chr") %>% 
  dplyr::filter(abs(distanceToTSS) < 3000) %>% 
  dplyr::select(seqnames, start, end, gene_symbol = SYMBOL) %>% 
  mutate(region = paste(seqnames, start, end, sep = "-"))
markers_for_cl0 = markers_for_cl0 %>% mutate(region = paste(chr, start, end, sep = "-")) %>% 
  inner_join(., markers_for_cl0_annot, by = "region") %>% 
  dplyr::select(region, gene_symbol, avg_log2FC, p_val, p_val_adj) %>% distinct_all()
markers_for_cl0 = markers_for_cl0 %>% dplyr::filter(p_val_adj < 0.05, avg_log2FC > 1) %>% mutate(cluster = 0)

# find G4 marker regions for cluster 1 (other brain cells)
markers_for_cl1 = FindMarkers(g4, test.use = "LR", latent.vars = "peak_region_fragments", ident.1 = 1, ident.2 = 0, group.by = "seurat_clusters") 
markers_for_cl1 = markers_for_cl1 %>% mutate(id = rownames(markers_for_cl1))
markers_for_cl1 = markers_for_cl1 %>% separate(id, sep = "-", into = c("chr", "start", "end"), remove = TRUE) %>% 
  mutate(start = as.numeric(start), end = as.numeric(end))
markers_for_cl1 = markers_for_cl1 %>% mutate(region = paste(chr, start, end, sep = "-"))
markers_for_cl1_annot = mm10_annotation(regions = markers_for_cl1, start_col = "start", end_col = "end", seqname_col = "chr") %>% 
  dplyr::filter(abs(distanceToTSS) < 3000) %>% 
  dplyr::select(seqnames, start, end, gene_symbol = SYMBOL) %>% 
  mutate(region = paste(seqnames, start, end, sep = "-"))
markers_for_cl1 = markers_for_cl1 %>% mutate(region = paste(chr, start, end, sep = "-")) %>% 
  inner_join(., markers_for_cl1_annot, by = "region") %>% 
  dplyr::select(region, gene_symbol, avg_log2FC, p_val, p_val_adj) %>% distinct_all()
markers_for_cl1 = markers_for_cl1 %>% dplyr::filter(p_val_adj < 0.05, avg_log2FC > 1) %>% mutate(cluster = 1)

# save
write_tsv(markers_for_cl1, glue("{result_folder}unsorted_brain_res0.1-G4markers_logreg-cluster1.tsv"))
write_tsv(markers_for_cl0, glue("{result_folder}unsorted_brain_res0.1-G4markers_logreg-cluster0.tsv"))

# making volcano
markers_for_cl1 = read_tsv(glue("{result_folder}unsorted_brain_res0.1-G4markers_logreg-cluster1.tsv"))
markers_for_cl0 = read_tsv(glue("{result_folder}unsorted_brain_res0.1-G4markers_logreg-cluster0.tsv"))

# add color, size, alpha indications
cols = c("up" = "#bdbdbd", "down" = "#addd8e", "unaltered" = "white")
sizes = c("up" = 4, "down" = 4, "unaltered" = 2)
alphas = c("up" = 1, "down" = 1, "unaltered" = 0.5)


markers_for_cl0 = markers_for_cl0 %>% mutate(avg_log2FC = avg_log2FC * (-1))
volc_input_res0.1 = rbind(markers_for_cl0, markers_for_cl1)
volc_input_res0.1 = volc_input_res0.1 %>% 
  mutate(group = case_when(
    avg_log2FC > 2 & p_val_adj < 0.05 ~ "up",
    avg_log2FC < -2 & p_val_adj < 0.05 ~ "down",
    avg_log2FC >= -2 & avg_log2FC <= 2 ~ "unaltered",
    TRUE ~ "non sign.")
  )
volc_input_res0.1 = volc_input_res0.1 %>% mutate(sign_label = case_when(
  avg_log2FC > 2 & p_val_adj < 1e-7 ~ gene_symbol,
  avg_log2FC < -2 & p_val_adj < 1e-7 ~ gene_symbol,
  avg_log2FC < -8 | avg_log2FC > 8 ~ gene_symbol,
  avg_log2FC >= -2 & avg_log2FC <= 2 ~ "",
  TRUE ~ ""
)) %>% mutate(cluster = as.character(cluster))

labels = volc_input_res0.1 %>% pull(sign_label) 

# plot
ggplot_volc_res0.1_fma = volc_input_res0.1 %>%
  ggplot(aes(x = avg_log2FC,
             y = -log10(p_val_adj),
             size = group,
             alpha = group,
             fill = group)) +
  geom_point(shape = 21, colour = "black") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = c(-2, 2),
             linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_shape_manual(values = c(21, 24)) +
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-10, 10, 2)),  	 
                     limits = c(-10, 10)) +
  labs(
    title = "Differential G-quadruplexed promoters",
    subtitle = "cluster 1 vs. cluster 0",
    x = "log2FoldChange",
    y = "-log10 adj.p-value",
    fill = "cluster",
    shape = ""
  ) +
  ylim(0, 30) +
  guides(alpha = FALSE, size = FALSE, fill = guide_legend(override.aes = list(size = 5))) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) + geom_text_repel(label = labels, size = 6, max.overlaps = 50) # add labels
ggplot_volc_res0.1_fma

ggsave(
  glue("{result_folder}unsorted_brain_res0.1-G4markers_logreg-volc.png"),
  plot = last_plot(),
  width = 9,
  height = 7,
  dpi = 300,
)

ggsave(
  glue("{result_folder}unsorted_brain_res0.1-G4markers_logreg-volc.pdf"),
  device = "pdf",
  plot = last_plot(),
  width = 9,
  height = 7,
  dpi = 300,
)

