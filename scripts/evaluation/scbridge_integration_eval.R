if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(
  "Seurat",
  "glue",
  "tidyverse",
  "data.table",
  "ggplot2",
  "ggpubr",
  "anndata",
  "zellkonverter",
  "RColorBrewer",
  "viridis",
  "wesanderson",
  "treemapify",
  "ComplexHeatmap",
  "circlize"
)

# result folder
result_folder = "../../results/scBridge/output/"

# helper function for cluster overlap analysis
source("../../utils/overlap_score.R")

# Seurat integrated object
g4 = readRDS(
  "../../results/Seurat/GFP_sorted_mousebrain/res0.8/integration/outputs/G4_scRNA_integration.Rds"
)

# read anndata objects (scBridge outputs)
# Bartosovic - GFP sorted G4 scCUT&Tag integration
scbr_g4 = readH5AD("../../results/scBridge/output/scCutTag_GFPsorted_gene_activity_scores-integrated.h5ad")
scbr_g4 = as.Seurat(scbr_g4, counts = "X", data = NULL)
scbr_rna = readH5AD("../../results/scBridge/output/scRNA_Seq-Bartosovic-integrated.h5ad")
scbr_rna = as.Seurat(scbr_rna, counts = "X", data = NULL)
comb = readH5AD("../../results/scBridge/output/combined-Bartosovic_GFPsorted.h5ad")
comb = as.Seurat(comb, counts = "X", data = NULL)

# Zeisel - unsorted (cluster1) G4 scCUT&Tag integration
scbr_g4_unsorted_cl1 = readH5AD(
  "../../results/scBridge/output/scCutTag_unsorted_cl1_gene_activity_scores-integrated.h5ad"
)
scbr_g4_unsorted_cl1 = as.Seurat(scbr_g4_unsorted_cl1, counts = "X", data = NULL)
scbr_rna_zeisel = readH5AD("../../results/scBridge/output/scRNA_Seq-Zeisel_et_al-neuron-integrated.h5ad")
scbr_rna_zeisel = as.Seurat(scbr_rna_zeisel, counts = "X", data = NULL)
comb_zeisel_uns_cl1 = readH5AD("../../results/scBridge/output/combined-Zeisel_neuron_unsorted_cl1.h5ad")
comb_zeisel_uns_cl1 = as.Seurat(comb_zeisel_uns_cl1, counts = "X", data = NULL)

# scBridge output tables
rel = fread("../../results/scBridge/output/Bartosovic_GFPsorted-scbridge_reliability.csv")
pred = fread("../../results/scBridge/output/Bartosovic_GFPsorted-scbridge_predictions.csv", header = TRUE)

scbr_rna@meta.data = scbr_rna@meta.data %>%
  mutate(CellType = str_replace_all(CellType, pattern = "Astrocytes", replacement = "AST")) %>%
  mutate(CellType = str_replace_all(CellType, pattern = "Oligodendrocytes", replacement = "MOL"))

scbr_g4@meta.data = scbr_g4@meta.data %>%
  mutate(Prediction = str_replace_all(Prediction, pattern = "Astrocytes", replacement = "AST")) %>%
  mutate(Prediction = str_replace_all(Prediction, pattern = "Oligodendrocytes", replacement = "MOL"))

glue(
  "number of AST: {as.character(scbr_g4@meta.data %>% dplyr::filter(Prediction == 'AST') %>% rownames %>% length)}
     number of non-AST: {as.character(scbr_g4@meta.data %>% dplyr::filter(Prediction != 'AST') %>% rownames %>% length)}"
)

glue(
  "Av. reliability of AST: {as.character(scbr_g4@meta.data %>% dplyr::filter(Prediction == 'AST') %>% pull('Reliability') %>%
  mean %>% round(2))}"
)

glue(
  "# of reliable cells: {as.character(scbr_g4@meta.data %>% dplyr::filter(Reliability > 0.9) %>% rownames %>%
  length)}
  # of all cells: {as.character(dim(scbr_g4@meta.data)[1])}"
)

### scBridge integration of sorted G4 scCUT&Tag and Bartosovic scRNA-Seq ###
# Seurat UMAPs
set3 = brewer.pal(8, "Set3")
cols = c(
  'AST' = set3[1],
  'COP-NFOL' = set3[2],
  'MOL' = set3[3],
  'OPC' = set3[4],
  'OEC' = set3[5],
  'VEC' = set3[6],
  'VLMC' = set3[7],
  'Pericytes' = set3[8]
)

# cell types
scbr_rna_pred = DimPlot(
  object = scbr_rna,
  group.by = "CellType",
  label = FALSE,
  pt.size = 2,
  label.size = 7,
  raster = FALSE,
  repel = TRUE,
  cols = cols,
  order = c(
    'Pericytes',
    'VLMC',
    'VEC',
    'COP-NFOL',
    'OPC',
    'OEC',
    'MOL',
    'AST'
  )
) +
  xlim(-12, 25) +
  ylim(-20, 20) +
  ggtitle("Cell type (scRNA-Seq)") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )
scbr_rna_pred

scbr_g4@meta.data = scbr_g4@meta.data %>%
  mutate(Prediction = as.character(Prediction)) %>%
  mutate(Prediction = ifelse(str_detect(Prediction, "Novel"), "unreliable", Prediction)) %>%
  mutate(Prediction = as.factor(Prediction))

cols = c(
  'AST' = set3[1],
  'COP-NFOL' = set3[2],
  'MOL' = set3[3],
  'OPC' = set3[4],
  'OEC' = set3[5],
  'VEC' = set3[6],
  'VLMC' = set3[7],
  'Pericytes' = set3[8],
  'unreliable' = '#f0f0f0'
)

# predicted labels by scBridge
scbr_g4_pred = DimPlot(
  object = scbr_g4,
  group.by = "Prediction",
  label = FALSE,
  pt.size = 2,
  label.size = 7,
  repel = TRUE,
  raster = FALSE,
  cols = cols,
  order = c(
    'Pericytes',
    'VLMC',
    'VEC',
    'COP-NFOL',
    'OPC',
    'OEC',
    'MOL',
    'AST',
    'unreliable'
  )
) +
  xlim(-12, 25) +
  ylim(-20, 20) +
  ggtitle("Prediction") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )
scbr_g4_pred

# scBridge reliability
scbr_rel = FeaturePlot(object = scbr_g4,
                  features = 'Reliability',
                  raster = FALSE) +
  scale_color_viridis() +
  xlim(-12, 25) +
  ylim(-20, 20) +
  ggtitle("Reliability") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )
scbr_rel

# visualize Seurat prediction scores
g4@meta.data = g4@meta.data %>% rownames_to_column(., var = "cell_id")
scbr_g4@meta.data = scbr_g4@meta.data %>% rownames_to_column(., var = "cell_id")
scbr_g4@meta.data = scbr_g4@meta.data %>% left_join(., g4@meta.data, by = "cell_id")
rownames(scbr_g4@meta.data) = scbr_g4@meta.data$cell_id
scbr_g4 = subset(scbr_g4, cells = g4@meta.data$cell_id)

scbr_g4@meta.data = scbr_g4@meta.data %>%
  dplyr::select(
    cell_id,
    Domain,
    Prediction,
    Reliability,
    Seurat_prediction = pred_cell_type,
    Seurat_pred_score = pred_max_score
  )

seurat_pred_score = FeaturePlot(object = scbr_g4,
                                features = 'Seurat_pred_score',
                                raster = FALSE) +
  scale_color_viridis() +
  xlim(-12, 25) +
  ylim(-20, 20) +
  ggtitle("Seurat prediction score") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )
seurat_pred_score

scbr_g4@meta.data = scbr_g4@meta.data %>%
  mutate(Seurat_prediction =
           str_replace_all(
             Seurat_prediction,
             pattern = "Astrocytes",
             replacement = "AST"
           )) %>%
  mutate(
    Seurat_prediction =
      str_replace_all(
        Seurat_prediction,
        pattern = "Oligodendrocytes",
        replacement = "MOL"
      )
  )

cols = c(
  'AST' = set3[1],
  'COP-NFOL' = set3[2],
  'MOL' = set3[3],
  'OPC' = set3[4],
  'OEC' = set3[5],
  'VEC' = set3[6],
  'VLMC' = set3[7],
  'Pericytes' = set3[8]
)

seurat_pred = DimPlot(
  object = scbr_g4,
  group.by = "Seurat_prediction",
  label = FALSE,
  pt.size = 2,
  label.size = 7,
  repel = TRUE,
  raster = FALSE,
  cols = cols,
  order = c(
    'Pericytes',
    'VLMC',
    'VEC',
    'COP-NFOL',
    'OPC',
    'OEC',
    'MOL',
    'AST'
  )
) +
  xlim(-12, 25) +
  ylim(-20, 20) +
  ggtitle("Seurat prediction") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )
seurat_pred

scbr_g4@meta.data = scbr_g4@meta.data %>% mutate(AST_status = ifelse(Prediction == "AST", "AST", "non-AST"))
set3 = brewer.pal(9, "Set3")
cols = c("AST" = set3[1], "non-AST" = set3[9])

ast = DimPlot(
  object = scbr_g4,
  group.by = "AST_status",
  label = FALSE,
  pt.size = 2,
  label.size = 7,
  repel = TRUE,
  raster = FALSE,
  order = c("AST", "non-AST")
) +
  xlim(-12, 25) +
  ylim(-20, 20) +
  scale_color_manual(values = cols) +
  ggtitle("") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )
ast

ggsave(
  glue("{result_folder}Seurat_predAST-UMAP.pdf"),
  plot = ast,
  width = 8,
  height = 8,
  device = "pdf"
)

# experimental system (modality)
comb@meta.data = comb@meta.data %>%
  mutate(Domain = as.character(Domain)) %>%
  mutate(Domain = ifelse(str_detect(Domain, "scCutTag_"), "G4 scCut&Tag", Domain)) %>%
  mutate(Domain = ifelse(str_detect(Domain, "Bartosovic"), "scRNA-Seq", Domain)) %>%
  mutate(Domain = as.character(Domain))

domain = DimPlot(
  object = comb,
  group.by = "Domain",
  label = FALSE,
  pt.size = 2,
  label.size = 7,
  repel = TRUE,
  raster = FALSE
) +
  xlim(-12, 25) +
  ylim(-20, 20) +
  ggtitle("Modality") +
  scale_fill_manual(values = c("#fc9272", "#a6bddb")) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )
domain

umaps = ggarrange(
  plotlist = list(
    scbr_rna_pred,
    scbr_g4_pred,
    scbr_rel,
    domain,
    seurat_pred_score,
    seurat_pred
  ),
  ncol = 2,
  nrow = 3
)
ggsave(
  glue("{result_folder}Seurat_UMAPs.pdf"),
  plot = umaps,
  width = 18,
  height = 18,
  device = "pdf"
)

featureplots = list()
genes = c("Tnik", "Pitpnc1", "Pbx1", "Nwd1")
for (gene in genes) {
  print(gene)
  plot = FeaturePlot(
    object = scbr_g4,
    features = gene,
    order = TRUE,
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
  featureplots[[gene]] = plot
}

ast_featureplotes = ggarrange(plotlist = featureplots,
                              ncol = 2,
                              nrow = 2)

ggsave(
  glue("{result_folder}Seurat_UMAPs-AST_spec_G4_markers.pdf"),
  plot = ast_featureplotes,
  width = 18,
  height = 18,
  device = "pdf"
)

featureplots_rna = list()
genes = c("Tnik", "Pitpnc1", "Pbx1", "Nwd1")
for (gene in genes) {
  plot = FeaturePlot(
    object = scbr_rna,
    features = gene,
    order = FALSE,
    raster = FALSE,
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

ggsave(
  glue(
    "{result_folder}Seurat_UMAPs-AST_spec_G4_markers-RNA_level.pdf"
  ),
  plot = ast_featureplotes_rna,
  width = 18,
  height = 18,
  device = "pdf"
)

### scBridge integration of unsorted G4 scCUT&Tag (cluster 1) and Zeisel neuron scRNA-Seq ###
glue(
  "number of reliable preds: \n {as.character(scbr_g4_unsorted_cl1@meta.data %>% 
  dplyr::filter(Reliability > 0.90) %>% rownames %>% length)}/{as.character(dim(scbr_g4_unsorted_cl1@meta.data)[1])}"
)

glue(
  "Av. reliability: {as.character(scbr_g4_unsorted_cl1@meta.data %>% pull('Reliability') %>%
  mean %>% round(2))}"
)

glue(
  "Number of the unreliable predictions:
  {as.character(scbr_g4_unsorted_cl1@meta.data %>% 
  dplyr::filter(str_detect(Prediction, 'Unreliable')) %>%
  rownames %>% length)}"
)


# cell types
cols2 = c(
  'amygdala' = "#fc9272",
  'cerebellum' = '#3182bd',
  'cortex1' = "#a1d99b",
  'cortex2' = "#31a354",
  'cortex3' = "#006d2c",
  "drg" = "#8c96c6",
  "enteric" = "#feb24c",
  "hippocampus" = "#ffffb2",
  "hypothalamus" = "#8c2d04",
  "medulla" = "#deebf7",
  "midbraindorsal" = "#e5f5e0",
  "midbrainventral" = "#636363",
  "olfactory" = "#016450",
  "pons" = "#c51b8a",
  "striatumdorsal" = "#e31a1c",
  "striatumventral" = "#9ecae1",
  "sympathetic" = "#8856a7",
  "thalamus" = "#810f7c"
)

scbr_rna_zeisel_pred = DimPlot(
  object = scbr_rna_zeisel,
  group.by = "CellType",
  label = FALSE,
  pt.size = 0.5,
  label.size = 7,
  raster = FALSE,
  repel = TRUE,
  cols = cols2
) +
  labs(color = "neuron type") +
  xlim(-10, 20) +
  ylim(-10, 20) +
  ggtitle("Cell type (scRNA-Seq)") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme(
    text = element_text(size = 10),
    plot.title = element_text(size = 20),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )
scbr_rna_zeisel_pred

# predicted labels by scBridge
cols3 = c(
  'Novel (Most Unreliable)' = '#f0f0f0',
  'amygdala' = "#fc9272",
  'cerebellum' = '#3182bd',
  'cortex1' = "#a1d99b",
  'cortex2' = "#31a354",
  'cortex3' = "#006d2c",
  "drg" = "#8c96c6",
  "enteric" = "#feb24c",
  "hippocampus" = "#ffffb2",
  "hypothalamus" = "#8c2d04",
  "medulla" = "#deebf7",
  "midbraindorsal" = "#e5f5e0",
  "midbrainventral" = "#636363",
  "olfactory" = "#016450",
  "pons" = "#c51b8a",
  "striatumdorsal" = "#e31a1c",
  "striatumventral" = "#9ecae1",
  "sympathetic" = "#8856a7",
  "thalamus" = "#810f7c"
)

scbr_g4_unsorted_cl1_pred = DimPlot(
  object = scbr_g4_unsorted_cl1,
  group.by = "Prediction",
  label = FALSE,
  pt.size = 0.5,
  label.size = 7,
  repel = TRUE,
  raster = FALSE,
  cols = cols3
) +
  labs(color = "pred. neuron type") +
  xlim(-10, 20) +
  ylim(-10, 20) +
  ggtitle("Prediction") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme(
    text = element_text(size = 10),
    plot.title = element_text(size = 20),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )
scbr_g4_unsorted_cl1_pred

# scBridge reliability
scbr_g4_unsorted_cl1_rel = FeaturePlot(object = scbr_g4_unsorted_cl1,
                                       features = 'Reliability',
                                       raster = FALSE) +
  scale_color_viridis() +
  xlim(-10, 20) +
  ylim(-10, 20) +
  ggtitle("Reliability") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme(
    text = element_text(size = 10),
    plot.title = element_text(size = 20),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )
scbr_g4_unsorted_cl1_rel

# experimental system (modality)
comb_zeisel_uns_cl1@meta.data = comb_zeisel_uns_cl1@meta.data %>%
  mutate(Domain = as.character(Domain)) %>%
  mutate(Domain = ifelse(
    str_detect(Domain, "G4 scCut&Tag"),
    "G4 scCut&Tag (unsorted, cl 1)",
    Domain
  )) %>%
  mutate(Domain = ifelse(
    str_detect(Domain, "scRNA_Seq"),
    "neuron scRNA_Seq (Zeisel et al.)",
    Domain
  )) %>%
  mutate(Domain = as.character(Domain))

domain_unsorted_cl1 = DimPlot(
  object = comb_zeisel_uns_cl1,
  group.by = "Domain",
  label = FALSE,
  pt.size = 0.5,
  label.size = 7,
  repel = TRUE,
  raster = FALSE,
  alpha = 1
) +
  xlim(-10, 20) +
  ylim(-10, 20) +
  ggtitle("Modality") +
  scale_fill_manual(values = c("#fc9272", "#a6bddb")) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme(
    text = element_text(size = 10),
    plot.title = element_text(size = 20),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )
domain_unsorted_cl1

umaps_unsorted_int = ggarrange(
  plotlist = list(
    scbr_rna_zeisel_pred,
    scbr_g4_unsorted_cl1_pred,
    domain_unsorted_cl1,
    scbr_g4_unsorted_cl1_rel
  ),
  ncol = 2,
  nrow = 2
)
umaps_unsorted_int
ggsave(
  glue("{result_folder}Seurat_UMAPs-unsorted_cl1_int.pdf"),
  plot = umaps_unsorted_int,
  width = 18,
  height = 18,
  device = "pdf"
)

meta = scbr_g4_unsorted_cl1@meta.data
unsorted_int_preds = meta %>%
  mutate(
    Prediction = case_when(
      Prediction == "Novel (Most Unreliable)" ~ "unreliable",
      str_detect(Prediction, "cortex") ~ "cortex",
      Prediction == "drg" ~ "dorsal root ganglia",
      TRUE ~ as.character(Prediction)
    )
  ) %>%
  group_by(Prediction) %>% dplyr::count() %>%
  dplyr::rename("count" = n) %>%
  mutate(log_count = log2(count))

# visualizing the predicted cell types of unsorted cluster 1
ggplot(unsorted_int_preds,
       aes(area = log_count, fill = log_count, label = Prediction)) +
  geom_treemap() +
  geom_treemap_text(colour = "white",
                    place = "centre",
                    size = 15) +
  labs(fill = "log2 count") +
  scale_fill_viridis_c()

ggsave(
  glue("{result_folder}tree_plot-unsorted_cl1_predictions.pdf"),
  plot = last_plot(),
  width = 6,
  height = 6,
  device = "pdf"
)

# cluster overlaps (overlap score heatmaps)
# scBridge clusters - Seurat clusters
# scBridge integration
col_fun = colorRamp2(c(0, 0.25, 0.5), c("#452258", "#679b81", "#f0e527"))

meta = g4@meta.data %>%
  left_join(., rel, c("cell_id" = "V1")) %>% dplyr::rename(scBridge_reliability = "Reliability") %>%
  left_join(., pred, c("cell_id" = "V1")) %>% dplyr::rename(scBridge_prediction = "Prediction") %>%
  filter(scBridge_reliability > 0.9) %>% dplyr::filter(scBridge_prediction != "Novel (Most Unreliable)") %>%
  dplyr::filter(pred_max_score > 0.5) %>%
  mutate(
    scBridge_prediction = str_replace_all(
      scBridge_prediction,
      pattern = "Astrocytes",
      replacement = "AST"
    )
  ) %>%
  mutate(
    scBridge_prediction = str_replace_all(
      scBridge_prediction,
      pattern = "Oligodendrocytes",
      replacement = "MOL"
    )
  ) %>%
  mutate(pred_cell_type = str_replace_all(pred_cell_type, pattern = "Astrocytes", replacement = "AST")) %>%
  mutate(pred_cell_type = str_replace_all(pred_cell_type, pattern = "Oligodendrocytes", replacement = "MOL"))

ident2seurat = data.frame(idents = rownames(meta), rna_label = meta$pred_cell_type)
ident2seurat = ident2seurat[complete.cases(ident2seurat), ]

ident2sc_br = data.frame(idents = rownames(meta), rna_label = meta$scBridge_prediction)
ident2sc_br = ident2sc_br[complete.cases(ident2sc_br), ]

ovlpScore.df = cal_ovlpScore(ident2seurat, ident2sc_br)

mapSubclass = ovlpScore.df
colnames(mapSubclass) <- c("cell_type", "seurat_cluster", "ovlpScore")
mapSubclass = dcast(mapSubclass, cell_type~seurat_cluster, value.var = "ovlpScore", 
                    fun.aggregate = identity, fill = 0)
rows = mapSubclass$cell_type
mapSubclass = mapSubclass[,-1]
rownames(mapSubclass) = rows

pdf(
  file = glue("{result_folder}scBridge_Seurat_overlap_score_hm.pdf"),
  width = 6,
  height = 5
)
overlap_hm_scbr_seurat = Heatmap(
  mapSubclass,
  name = "overlap score",
  clustering_distance_rows = "pearson",
  col = col_fun,
  row_title = "Seurat prediction (pred. score > 0.5)",
  row_title_side = "right",
  column_title = "scBridge prediction",
  column_title_side = "bottom",
  show_row_names = TRUE,
  rect_gp = gpar(col = "black", lwd = 1),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  width = unit(7, "cm"),
  height = unit(3, "cm"),
  column_names_rot = 90
)
overlap_hm_scbr_seurat
dev.off()

# scBridge clusters
# scBridge integration
ident2rna = data.frame(idents = rownames(meta), rna_label = meta$scBridge_prediction)
ident2rna = ident2rna[complete.cases(ident2rna), ]
ident2rna = ident2rna[which(ident2rna$rna_label != "Novel (Most Unreliable)"), ]

ident2g4 = data.frame(idents = rownames(meta), rna_label = as.character(meta$seurat_clusters))
ident2g4 = ident2g4[complete.cases(ident2g4), ]

ovlpScore.df = cal_ovlpScore(ident2rna, ident2g4)

mapSubclass = ovlpScore.df
colnames(mapSubclass) <- c("cell_type", "seurat_cluster", "ovlpScore")
mapSubclass = dcast(mapSubclass, cell_type~seurat_cluster, value.var = "ovlpScore", 
                    fun.aggregate = identity, fill = 0)
rows = mapSubclass$cell_type
mapSubclass = mapSubclass[,-1]
rownames(mapSubclass) = rows

pdf(
  file = glue("{result_folder}scBridge_overlap_score_hm.pdf"),
  width = 6,
  height = 4
)
overlap_hm_scbr = Heatmap(
  mapSubclass,
  name = "overlap score",
  clustering_distance_rows = "pearson",
  col = col_fun,
  row_title = "scBridge prediction",
  row_title_side = "right",
  column_title = "Seurat cluster",
  column_title_side = "bottom",
  show_row_names = TRUE,
  rect_gp = gpar(col = "black", lwd = 1),
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  width = unit(7, "cm"),
  height = unit(7, "cm"),
  column_names_rot = 0
)
overlap_hm_scbr
dev.off()
