if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(
  "tidyverse",
  "Seurat",
  "zellkonverter",
  "Signac",
  "wigglescout",
  "data.table",
  "ggplot2",
  "ggpubr",
  "glue",
  "RColorBrewer",
  "ggrepel",
  "ggrastr",
  "viridis"
)

# helper function for annotation
source("../../utils/annotation.R")

# folders
result_folder = "../../results/astrocyte_exploration/"

# Seurat objects
sorted = readRDS(file = "../../results/Seurat/GFP_sorted_mousebrain/res0.8/outputs/Seurat_object.Rds")
rna = readRDS(file = "../../data/scRNA-Seq/scRNA_Seq-mouse_brain.Rds")
rna@meta.data = rna@meta.data %>%
  mutate(cell_type = str_replace_all(cell_type, pattern = "Astrocytes", replacement = "AST")) %>%
  mutate(cell_type = str_replace_all(cell_type, pattern = "Oligodendrocytes", replacement = "MOL"))
norm = rna[['RNA']]@data

# scRNA-Seq markers
markers = read_tsv(
  "../../results/Seurat/GFP_sorted_mousebrain/res0.8/integration/outputs/scRNA-Seq-FindAllMarkers_output.tsv"
)
markers = markers[markers$p_val < 0.05 &
                    markers$avg_log2FC > 0.5, ]
markers = markers %>%
  mutate(cluster = str_replace_all(cluster, pattern = "Astrocytes", replacement = "AST")) %>%
  mutate(cluster = str_replace_all(cluster, pattern = "Oligodendrocytes", replacement = "MOL"))

# scBridge outputs
scbr_rna = readH5AD("../../results/scBridge/output/Bartosovic_scRNA-Seq-integrated.h5ad")
scbr_rna = as.Seurat(scbr_rna, counts = "X", data = NULL)

rel = fread("../../results/scBridge/output/scbridge_reliability.csv")
pred = fread("../../results/scBridge/output/scbridge_predictions.csv",
             header = TRUE)
pred = pred %>%
  mutate(Prediction = str_replace_all(Prediction, pattern = "Astrocytes", replacement = "AST")) %>%
  mutate(Prediction = str_replace_all(Prediction, pattern = "Oligodendrocytes", replacement = "MOL"))

# Seurat prediction scores
pred_score = readRDS(
  "../../results/Seurat/GFP_sorted_mousebrain/res0.8/integration/outputs/g4_cell_label_preds.Rds"
)

pred_score = t(pred_score@data)
ids = rownames(pred_score)
pred_score = as_tibble(pred_score)
colnames(pred_score) = c("OEC",
                         "AST",
                         "MOL",
                         "Pericytes",
                         "VEC",
                         "VLMC",
                         "COP-NFOL",
                         "OPC",
                         "max")
pred_score = pred_score %>%
  dplyr::select(-max) %>% mutate(cell_id = ids) %>% dplyr::filter(AST > 0.90)
ast_ids = pred_score %>% pull(cell_id)

# AST scRNA-Seq markers
ast_markers = markers %>%
  mutate(cluster = ifelse(str_detect(markers$cluster, "AST"), "AST", cluster)) %>%
  dplyr::filter(str_detect(cluster, "AST")) %>%
  arrange(desc(avg_log2FC)) %>%
  top_n(500, wt = avg_log2FC)
all_ast_markers = ast_markers$gene
ast_markers = unique(ast_markers$gene)[1:20]

# highly predicted AST cell ids
barcodes_seurat = pred_score$cell_id[which(pred_score$AST > 0.90)]
barcodes_scbridge = pred %>% filter(Prediction == "AST") %>% pull(V1)
barcodes_scbridge_others = pred %>% filter(Prediction != "AST") %>%
  filter(Prediction != "Novel (Most Unreliable)") %>%
  pull(V1)

barcodes = unique(c(barcodes_seurat, barcodes_scbridge))
barcodes_df = tibble(AST = barcodes, type = "AST")
write_tsv(barcodes_df, glue("{result_folder}Pred_AST-barcodes.tsv"))

barcodes_non_ast = rownames(sorted@meta.data)[which(!rownames(sorted@meta.data) %in% barcodes)]
barcodes_non_ast_df = tibble(non_AST = barcodes_non_ast, type = "non_AST")
write_tsv(barcodes_non_ast_df,
          glue("{result_folder}Pred_non_AST-barcodes.tsv"))

# differential GA scores
barcodes = barcodes_df %>% dplyr::filter(AST %in% rownames(sorted@meta.data)) %>%
  pull(AST)
barcodes_non_ast = barcodes_non_ast_df %>% dplyr::filter(non_AST %in% rownames(sorted@meta.data)) %>%
  pull(non_AST)

diff_AST_vs_nonAST_seurat = FindMarkers(
  sorted,
  ident.1 = barcodes,
  ident.2 = barcodes_non_ast,
  only.pos = FALSE,
  assay = "GA",
  logfc.threshold = 0,
  test.use = "LR",
  latent.vars = "peak_region_fragments"
)

write.table(
  diff_AST_vs_nonAST_seurat,
  glue("{result_folder}diff_AST_vs_nonAST_seurat.tsv"),
  quote = FALSE
)

sign_AST_vs_nonAST_seurat = diff_AST_vs_nonAST_seurat %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  arrange(desc(avg_log2FC)) %>%
  rownames

barcodes_scbridge = pred %>%
  dplyr::filter(Prediction == "AST") %>%
  dplyr::filter(V1 %in% rownames(sorted@meta.data)) %>%
  pull(V1)
barcodes_scbridge_others = pred %>%
  dplyr::filter(Prediction != "AST") %>%
  dplyr::filter(V1 %in% rownames(sorted@meta.data)) %>%
  pull(V1)

diff_AST_vs_nonAST_scBr = FindMarkers(
  sorted,
  ident.1 = barcodes_scbridge,
  ident.2 = barcodes_scbridge_others,
  only.pos = FALSE,
  assay = "GA",
  logfc.threshold = 0,
  test.use = "LR",
  latent.vars = "peak_region_fragments"
)

write.table(
  diff_AST_vs_nonAST_scBr,
  glue("{result_folder}diff_AST_vs_nonAST_scBridge.tsv"),
  quote = FALSE
)

sign_AST_vs_nonAST_scBr = diff_AST_vs_nonAST_scBr %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  arrange(desc(avg_log2FC)) %>%
  rownames

## visualizations
meta = sorted@meta.data
meta = meta %>%
  mutate(AST_status = ifelse(
    rownames(meta) %in% barcodes,
    "predicted AST",
    "predicted non-AST"
  )) %>%
  mutate(
    AST_status_scBridge = ifelse(
      rownames(meta) %in% barcodes_scbridge,
      "predicted AST",
      "predicted non-AST"
    )
  )
sorted@meta.data = meta

genes = markers %>% dplyr::filter(cluster == "AST") %>% arrange(desc(avg_log2FC)) %>% top_n(6, wt = avg_log2FC) %>% pull(gene)
set3 = brewer.pal(9, "Set3")
cols = c("AST" = set3[1], "non-AST" = set3[9])

plots = list()
for (i in seq(1:length(genes))) {
  plots[[i]] = CoveragePlot(
    object = sorted,
    region = genes[i],
    annotation = TRUE,
    show.bulk = TRUE,
    ymax = 10,
    group.by = "AST_status",
    peaks = TRUE
  )
  print(plot)
}

examples = ggarrange(plotlist = plots,
                     ncol = 2,
                     nrow = 3)

ggsave(
  glue("{result_folder}scBr_predictions-browser_example.pdf"),
  plot = examples,
  width = 12,
  height = 12,
  device = "pdf"
)

ggsave(
  glue("{result_folder}scBr_predictions-browser_example.png"),
  plot = examples,
  width = 12,
  height = 12
)

plots = list()
for (i in seq(1:length(genes))) {
  mol = norm[genes[i], rownames(rna@meta.data[which(rna@meta.data$cell_type == "AST"), ])]
  mol = tibble(type = "AST", expr = mol)
  nonmol = norm[genes[i], rownames(rna@meta.data[which(rna@meta.data$cell_type != "non-AST"), ])]
  nonmol = tibble(type = "non-AST", expr = nonmol)
  both = rbind(mol, nonmol)
  plots[[i]] = ggplot(both, aes(x = type, y = expr, fill = type)) +
    geom_boxplot(color = "black") +
    ylim(0, 5) +
    labs(title = genes[i], x = "", y = "log norm. expr") +
    scale_fill_manual(values = cols) +
    theme_classic() +
    theme(
      text = element_text(size = 20),
      plot.title = element_text(size = 15),
      axis.text.x = element_text(size = 20, color = "black"),
      axis.text.y = element_text(size = 20, color = "black")
    ) +
    stat_compare_means(label.y = 160,
                       label.x = 1.15,
                       size = 6)
}

bp_examples = ggarrange(plotlist = plots,
                        ncol = 2,
                        nrow = 3)

ggsave(
  glue("{result_folder}AST_markers-browser_example-expr_bps.pdf"),
  plot = bp_examples,
  width = 12,
  height = 12,
  device = "pdf"
)

ggsave(
  glue("{result_folder}AST_markers-browser_example-expr_bps.png"),
  plot = bp_examples,
  width = 12,
  height = 12
)

plots = list()
for (i in seq(1:length(sign_AST_vs_nonAST_scBr[1:10]))) {
  ast = norm[sign_AST_vs_nonAST_scBr[i], rownames(rna@meta.data[which(rna@meta.data$cell_type == "AST"), ])]
  ast = tibble(type = "AST", expr = ast)
  nonast = norm[sign_AST_vs_nonAST_scBr[i], rownames(rna@meta.data[which(rna@meta.data$cell_type != "AST"), ])]
  nonast = tibble(type = "non-AST", expr = nonast)
  both = rbind(ast, nonast)
  plots[[i]] = ggplot(both, aes(x = type, y = expr, fill = type)) +
    geom_boxplot(color = "black") +
    ylim(0, 6) +
    labs(title = sign_AST_vs_nonAST_scBr[i], x = "", y = "log norm. expr") +
    scale_fill_manual(values = cols) +
    theme_classic() +
    theme(
      text = element_text(size = 20),
      plot.title = element_text(size = 15),
      axis.text.x = element_text(size = 20, color = "black"),
      axis.text.y = element_text(size = 20, color = "black")
    ) +
    stat_compare_means(label.y = 7,
                       label.x = 1.15,
                       size = 6)
}

sign_scBridge_genes_expr_bp = ggarrange(plotlist = plots,
                                        ncol = 2,
                                        nrow = 5)
sign_scBridge_genes_expr_bp

ggsave(
  glue(
    "{result_folder}scBr_AST_markers-browser_example-expr_bps.png"
  ),
  plot = sign_scBridge_genes_expr_bp,
  width = 12,
  height = 16
)

ggsave(
  glue(
    "{result_folder}scBr_AST_markers-browser_example-expr_bps.pdf"
  ),
  plot = sign_scBridge_genes_expr_bp,
  width = 12,
  height = 16,
  device = "pdf"
)

# volcano plot AST vs. non-AST
volc_input = diff_AST_vs_nonAST_scBr %>%
  mutate(gene_name = rownames(.))

volc_input = volc_input %>% mutate(
  group = case_when(
    avg_log2FC > 0.5 & p_val_adj < 0.05 ~ "up",
    avg_log2FC < -0.5 & p_val_adj < 0.05 ~ "down",
    avg_log2FC >= -0.5 & avg_log2FC <= 0.5 ~ "unaltered",
    p_val_adj > 0.05 ~ "unaltered"
  )
) %>%
  mutate(
    sign_label = case_when(
      avg_log2FC > 0.5 & p_val_adj < 0.05 ~ gene_name,
      avg_log2FC < -0.5 &
        p_val_adj < 0.05 ~ gene_name,
      avg_log2FC >= -0.5 &
        avg_log2FC <= 0.5 ~ ""
    )
  )
labels = volc_input %>% pull(sign_label)

cols = c(
  "up" = brewer.pal(9, "Set3")[1],
  "down" = brewer.pal(9, "Set3")[9],
  "unaltered" = "white"
)
sizes = c("up" = 4,
          "down" = 4,
          "unaltered" = 2)
alphas = c("up" = 1,
           "down" = 1,
           "unaltered" = 0.5)

glue(
  "sign. down: {as.character(
     volc_input %>% dplyr::filter(group == 'down') %>% rownames %>% length)}
     sign. up: {as.character(
     volc_input %>% dplyr::filter(group == 'up') %>% rownames %>% length)}"
)

ggplot_volc = volc_input %>%
  ggplot(aes(
    x = avg_log2FC,
    y = -log10(p_val_adj),
    fill = group,
    size = group,
    alpha = group
  )) +
  ggrastr::geom_point_rast(shape = 21, # Specify shape and colour as fixed local parameters
                           colour = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-6, 6, 2)), limits = c(-6, 6)) +
  scale_y_continuous(breaks = c(seq(0, 10, 2)), limits = c(0, 10)) +
  labs(title = "predicted AST vs other G4 differences",
       x = "log2FoldChange",
       y = "-log10 adj.p-value",
       fill = " ") +
  guides(alpha = FALSE,
         size = FALSE,
         fill = guide_legend(override.aes = list(size = 5))) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) +
  geom_text_repel(label = labels,
                  size = 7,
                  max.overlaps = 100) # add labels
ggplot_volc

ggsave(
  glue("{result_folder}scBr_predAST_vs_nonAST_volc.pdf"),
  width = 10,
  height = 10,
  device = "pdf"
)
ggsave(
  glue("{result_folder}scBr_predAST_vs_nonAST_volc.png"),
  width = 10,
  height = 10,
  dpi = 300
)

# feature plots of scBridge AST markers
featureplots_rna = list()
genes = c("Tnik", "Pitpnc1", "Pbx1", "Nwd1")
for (gene in genes) {
  plot = FeaturePlot(
    object = scbr_rna,
    features = gene,
    order = FALSE,
    raster = TRUE,
    reduction = "X_umap",
    pt.size = 2
  ) +
    #scale_color_gradient2(low = "#f0f0f0", mid = "#f0f0f0", high = "red") +
    scale_color_viridis() +
    xlim(-12, 25) +
    ylim(-20, 20) +
    ggtitle(gene) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    theme(
      text = element_text(size = 25),
      plot.title = element_text(size = 20),
      axis.text.x = element_text(size = 25, color = "black"),
      axis.text.y = element_text(size = 25, color = "black")
    )
  featureplots_rna[[gene]] = plot
}

ast_featureplotes_rna = ggarrange(plotlist = featureplots_rna,
                                  ncol = 2,
                                  nrow = 2)
ast_featureplotes_rna

# export in pdf and png
ggsave(
  glue(
    "{result_folder}scBr_AST_markers-browser_example-expr_featureplots.pdf"
  ),
  width = 12,
  height = 10,
  device = "pdf"
)
ggsave(
  glue(
    "{result_folder}scBr_AST_markers-browser_example-expr_featureplots.png"
  ),
  width = 12,
  height = 10,
  dpi = 300
)
