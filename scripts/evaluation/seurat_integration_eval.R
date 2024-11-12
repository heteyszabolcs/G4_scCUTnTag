print("Load R packages")
if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(
  "Seurat",
  "Signac",
  "glue",
  "ggplot2",
  "dplyr",
  "tidyverse",
  "cowplot",
  "ggpubr",
  "circlize",
  "ComplexHeatmap"
)

# working directory of the integration
workdir = "../../results/Seurat/GFP_sorted_mousebrain/res0.8/integration/"

# load scRNA-Seq and results of Seurat integration
g4 = readRDS(glue("{workdir}/outputs/G4_scRNA_integration.Rds"))
rna = readRDS(glue("{workdir}/outputs/scRNA_Seq_Seurat_object.Rds"))
genes.use = VariableFeatures(rna)
cell_class_pred = readRDS(glue("{workdir}/outputs/g4_cell_label_preds.Rds"))

# load scRNA-Seq marker analysis
markers = read_tsv(glue("{workdir}/outputs/scRNA-Seq-FindAllMarkers_output.tsv"))

# visualizing predictions
# UMAP
pred_umap = DimPlot(
  g4,
  pt.size = 2,
  label.size = 7,
  repel = TRUE,
  raster = TRUE,
  group.by = "pred_cell_type"
) +
  scale_color_brewer(palette = "Pastel1") +
  xlim(-15, 15) +
  ylim(-15, 15) +
  ggtitle(" ") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) + NoAxes()
pred_umap

ggsave(
  glue("{workdir}/plots/predicted_labels_UMAP.pdf"),
  plot = pred_umap,
  width = 13,
  height = 10,
  device = "pdf"
)

# prediction score heatmap
predictions = as.matrix(cell_class_pred@data)
predictions = predictions[rownames(predictions) != "max", ]
col_fun = colorRamp2(c(0, 0.5, 1), c("#421654", "#458f8a", "#f0e527"))
pdf(
  file = glue("{workdir}/plots/predicton_score_heatmap.pdf"),
  width = 8,
  height = 6
)
Heatmap(
  predictions,
  column_title = " ",
  row_title = " ",
  name = "pred. score",
  # row_km = 2,
  # column_km = 2,
  clustering_method_rows = "complete",
  col = col_fun,
  #rect_gp = gpar(col = "black", lwd = 0.1),
  show_column_dend = TRUE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  show_column_names = FALSE,
  heatmap_width = unit(17, "cm"),
  heatmap_height = unit(6, "cm"),
  row_names_gp = gpar(fontsize = 12),
  column_names_rot = 90
)
dev.off()

# prediction score boxplot
types = rownames(predictions)
predictions = as.data.frame(predictions)
predictions$types = types
predictions_long = predictions %>% pivot_longer(
  .,
  cols = c("AAACGAAAGAAGCCGT-1":"TTTGTGTTCTCGCGTT-1"),
  names_to = "cell_id",
  values_to = "pred_score"
)
meta = as_tibble(g4@meta.data)
meta = meta %>% mutate(cell_id = rownames(g4@meta.data))
meta = predictions_long %>% inner_join(., meta, by = "cell_id")

pred_boxplots = lapply(types, function(x) {
  meta = meta %>% dplyr::filter(types == x)
  plot = ggplot(meta,
                aes(x = seurat_clusters, y = pred_score, fill = seurat_clusters)) +
    geom_boxplot(color = "black") +
    scale_fill_brewer(palette = "Pastel1") +
    ylim(0, 1) +
    labs(
      title = x,
      x = "",
      y = "",
      fill = ""
    ) +
    theme_classic() +
    guides(fill = "none") +
    theme(
      text = element_text(size = 9),
      plot.title = element_text(size = 15, face = "bold"),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.text.y = element_text(size = 12, color = "black")
    ) + stat_compare_means(label = "p.signif", label.y = 0.9)
  return(print(plot))
})

pred_boxplots = ggarrange(plotlist = pred_boxplots)

ggsave(
  glue("{workdir}/plots/predicton_score_boxplots.pdf"),
  plot = pred_boxplots,
  width = 10,
  height = 7,
  device = "pdf"
)

# coembedding G4 and scRNA-Seq
rna@meta.data$data_type = "scRNA-Seq"
g4@meta.data$data_type = "G4 scCut&Tag"

coembed = merge(x = rna, y = g4)
coembed = ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed = RunPCA(coembed, features = genes.use, verbose = FALSE)
ElbowPlot(coembed)
coembed = RunUMAP(coembed, dims = 1:15)

coembed.g4 = coembed[, coembed$data_type == "G4 scCut&Tag"]
coembed.scrna = coembed[, coembed$data_type == "scRNA-Seq"]

# functions for feature plots
create_expr_feature_plot = function(marker_gene) {
  # feature plot
  plot = FeaturePlot(
    object = coembed.scrna,
    features = marker_gene,
    min.cutoff = min(coembed.scrna@assays$RNA@data[marker_gene, ]),
    max.cutoff = max(coembed.scrna@assays$RNA@data[marker_gene, ]),
    pt.size = 4,
    raster = TRUE,
    order = TRUE
  ) +
    xlim(-15, 15) +
    ylim(-10, 10) +
    scale_color_gradient2(
      low = "#edf8b1",
      mid = "#7fcdbb",
      high = "#225ea8",
      midpoint = mean(c(
        min(coembed.scrna@assays$RNA@data[marker_gene, ]),
        max(coembed.scrna@assays$RNA@data[marker_gene, ])
      ))
    ) +
    theme(
      legend.position = 'bottom',
      text = element_text(size = 15),
      plot.title = element_text(size = 20),
      axis.text.x = element_text(size = 25, color = "black"),
      axis.text.y = element_text(size = 25, color = "black")
    ) +
    NoAxes()
  return(print(plot))
}

create_g4_feature_plot = function(marker_gene) {
  plot = FeaturePlot(
    object = coembed.g4,
    features = marker_gene,
    min.cutoff = min(coembed.g4@assays$GA@data[marker_gene, ]),
    max.cutoff = max(coembed.g4@assays$GA@data[marker_gene, ]),
    pt.size = 4,
    raster = TRUE,
    order = TRUE
  ) +
    xlim(-15, 15) +
    ylim(-10, 10) +
    scale_color_gradient2(
      low = "#fee5d9",
      mid = "#fcae91",
      high = "#cb181d",
      midpoint = mean(c(
        min(coembed.g4@assays$GA@data[marker_gene, ]),
        max(coembed.g4@assays$GA@data[marker_gene, ])
      ))
    ) +
    theme(
      legend.position = 'bottom',
      text = element_text(size = 15),
      plot.title = element_text(size = 20),
      axis.text.x = element_text(size = 25, color = "black"),
      axis.text.y = element_text(size = 25, color = "black")
    ) +
    NoAxes()
  return(print(plot))
}

check_featureplot_warnings = function(gene) {
  if (!gene %in% rownames(coembed.g4@assays$GA@data))
  {
    return(FALSE)
  }
  tt <- tryCatch(
    create_g4_feature_plot(gene),
    error = function(e)
      e,
    warning = function(w)
      w
  )
  tt2 <- tryCatch(
    create_g4_feature_plot(gene),
    error = function(e)
      e,
    warning = function(w)
      w
  )
  
  if (is(tt, "warning")) {
    return(FALSE)
  } else {
    return(TRUE)
  }
  
}

get_best_marker = function(cell_type) {
  print(cell_type)
  marker = markers %>% dplyr::filter(str_detect(cluster, cell_type)) %>%
    arrange(avg_log2FC) %>% top_n(n = 1, wt = avg_log2FC) %>% pull(gene)
  
  if (marker %in% rownames(coembed.g4@assays$GA@data) &
      check_featureplot_warnings(marker)) {
    return(marker)
  } else {
    marker = markers %>% dplyr::filter(str_detect(cluster, cell_type)) %>%
      arrange(avg_log2FC) %>% top_n(n = 40, wt = avg_log2FC) %>% pull(gene)
    for (i in marker) {
      if (i %in% rownames(coembed.g4@assays$GA@data) &
          check_featureplot_warnings(i)) {
        return(i)
      }
    }
  }
}

get_underexpr_marker = function(cell_type) {
  print(cell_type)
  marker = markers %>% dplyr::filter(str_detect(cluster, cell_type)) %>%
    arrange(avg_log2FC) %>% top_n(n = -1, wt = avg_log2FC) %>% pull(gene)
  
  if (marker %in% rownames(coembed.g4@assays$GA@data) &
      check_featureplot_warnings(marker)) {
    return(marker)
  } else {
    marker = markers %>% dplyr::filter(str_detect(cluster, cell_type)) %>%
      arrange(avg_log2FC) %>% top_n(n = -40, wt = avg_log2FC) %>% pull(gene)
    for (i in marker) {
      if (i %in% rownames(coembed.g4@assays$GA@data) &
          check_featureplot_warnings(i)) {
        return(i)
      }
    }
  }
}

cell_types = unique(coembed@meta.data$cell_type)
cell_types = cell_types[!is.na(cell_types)]

underexpr_scrna_markers = sapply(cell_types, get_underexpr_marker)
best_scrna_markers = sapply(cell_types, get_best_marker)
g4_feature_plots = lapply(best_scrna_markers, create_g4_feature_plot)
g4_feature_plots = ggarrange(plotlist = g4_feature_plots)
expr_feature_plots = lapply(best_scrna_markers, create_g4_feature_plot)
expr_feature_plots = ggarrange(plotlist = expr_feature_plots)

ggsave(
  glue("{workdir}/plots/G4_Feature_plots_scRNA-Seq_markers.pdf"),
  plot = g4_feature_plots,
  width = 10,
  height = 7,
  device = "pdf"
)

ggsave(
  glue("{workdir}/plots/expr_Feature_plots_scRNA-Seq_markers.pdf"),
  plot = expr_feature_plots,
  width = 10,
  height = 7,
  device = "pdf"
)

g4_feature_plots = lapply(underexpr_scrna_markers, create_g4_feature_plot)
g4_feature_plots = ggarrange(plotlist = g4_feature_plots)
expr_feature_plots = lapply(underexpr_scrna_markers, create_g4_feature_plot)
expr_feature_plots = ggarrange(plotlist = expr_feature_plots)

ggsave(
  glue(
    "{workdir}/plots/G4_Feature_plots_scRNA-Seq_neg_markers.pdf"
  ),
  plot = g4_feature_plots,
  width = 10,
  height = 7,
  device = "pdf"
)

ggsave(
  glue(
    "{workdir}/plots/expr_Feature_plots_scRNA-Seq_neg_markers.pdf"
  ),
  plot = expr_feature_plots,
  width = 10,
  height = 7,
  device = "pdf"
)

# Violin plots
violins = lapply(best_scrna_markers, function(x) {
  print(x)
  plot = VlnPlot(
    object = rna,
    features = x,
    raster = FALSE,
    log = FALSE,
    split.by = "cell_type",
    group.by = "cell_type"
  ) +
    ylab("Norm. expression level") +
    ylim(0, 4) +
    scale_fill_continuous(guide = FALSE) +
    scale_fill_manual(values = rep("#a6bddb", 8)) +
    xlab(" ") +
    guides(legend = "none", fill = "none") +
    theme(
      text = element_text(size = 15),
      plot.title = element_text(size = 20),
      axis.text.x = element_text(
        color = "black",
        angle = 90,
        vjust = 0.5
      ),
      axis.text.y = element_text(size = 11, color = "black"),
      axis.title.y = element_text(size = 11, color = "black")
    )
  print(plot)
  return(plot)
})

expr_violins = ggarrange(plotlist = violins)
expr_violins

ggsave(
  glue("{workdir}/plots/expr_Violin_plots_scRNA-Seq_markers.pdf"),
  plot = expr_violins,
  width = 8,
  height = 11,
  device = "pdf"
)

# heatmap of expression markers
hm = DoHeatmap(
  rna,
  features = best_scrna_markers,
  group.by = "cell_type",
  raster = FALSE,
  group.colors = brewer.pal(8, "Pastel1"),
  angle = 90
) +
  scale_fill_gradientn(colors = c("#421654", "#458f8a", "#f0e527")) +
  theme(axis.text.y = element_text(size = 18, color = "black")) + NoLegend()

ggsave(
  glue("{workdir}/plots/expr_heatmap_scRNA-Seq_markers.pdf"),
  plot = hm,
  width = 11,
  height = 10
)

## UMAP dimension plots
g4_umap = DimPlot(
  g4,
  pt.size = 2,
  label.size = 7,
  repel = TRUE,
  raster = TRUE
) +
  scale_color_brewer(palette = "Pastel1") +
  xlim(-15, 15) +
  ylim(-15, 15) +
  ggtitle(" ") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )
rna_umap = DimPlot(
  rna,
  pt.size = 2,
  label.size = 7,
  group.by = 'cell_type',
  repel = TRUE,
  raster = TRUE
) +
  scale_color_brewer(palette = "Set3") +
  xlim(-15, 15) +
  ylim(-15, 15) +
  ggtitle(" ") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )

coembed_cells = DimPlot(
  coembed,
  pt.size = 2,
  label.size = 7,
  group.by = 'cell_type',
  repel = TRUE,
  na.value = "grey50",
  raster = TRUE
) +
  scale_color_brewer(palette = "Pastel1") +
  xlim(-15, 15) +
  ylim(-15, 15) +
  ggtitle(" ") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )

coembed_cells2 = DimPlot(
  coembed,
  pt.size = 2,
  label.size = 7,
  group.by = 'cell_type',
  repel = TRUE,
  label = TRUE,
  na.value = "grey50",
  raster = TRUE
) +
  scale_color_brewer(palette = "Pastel1") +
  xlim(-15, 15) +
  ylim(-15, 15) +
  ggtitle(" ") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoLegend() + NoAxes()

ggsave(
  glue("{workdir}/plots/Seurat_scRNA-Seq_celltypes.pdf"),
  plot = coembed_cells2,
  width = 6,
  height = 6,
  device = "pdf"
)

coembed_clusters = DimPlot(
  coembed,
  pt.size = 2,
  label.size = 7,
  group.by = 'seurat_clusters',
  repel = TRUE,
  raster = TRUE
) +
  scale_color_brewer(palette = "Pastel1") +
  xlim(-15, 15) +
  ylim(-15, 15) +
  ggtitle(" ") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )

coembed_experiments = DimPlot(
  coembed,
  pt.size = 2,
  label.size = 7,
  group.by = 'data_type',
  repel = TRUE,
  na.value = "grey50",
  raster = TRUE
) +
  scale_colour_manual(values = c("#fc9272", "#9ecae1")) +
  xlim(-15, 15) +
  ylim(-15, 15) +
  ggtitle(" ") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )

coembed_predcelltype = DimPlot(
  coembed,
  pt.size = 2,
  label.size = 7,
  group.by = 'pred_cell_type',
  repel = TRUE,
  raster = TRUE
) +
  scale_color_brewer(palette = "Pastel1") +
  xlim(-15, 15) +
  ylim(-15, 15) +
  ggtitle(" ") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )

cells = coembed@meta.data %>% dplyr::filter(data_type == "G4 scCut&Tag") %>% rownames
coembed_predscore = FeaturePlot(
  coembed,
  pt.size = 2,
  cells = cells,
  label.size = 7,
  features = 'pred_max_score',
  raster = TRUE
) +
  xlim(-15, 15) +
  ylim(-15, 15) +
  ggtitle(" ") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )

coembed@meta.data$seurat_clusters[is.na(coembed@meta.data$seurat_clusters)] = "scRNA-Seq"
coembed_ps = coembed_cells + coembed_experiments + coembed_predcelltype + coembed_predscore
coembed_ps

ggsave(
  glue("{workdir}/plots/Seurat_scRNA-Seq_coembeds.pdf"),
  plot = coembed_ps,
  width = 12,
  height = 8,
  device = "pdf"
)

ggsave(
  glue(
    "{workdir}/plots/Seurat_scRNA-Seq_integration_datatypes.pdf"
  ),
  plot = coembed_experiments,
  width = 6,
  height = 6,
  device = "pdf"
)

umaps = plot_grid(g4_umap,
                  rna_umap,
                  coembed_experiments,
                  coembed_cells,
                  ncol = 2,
                  nrow = 2)

ggsave(
  glue("{workdir}/plots/Seurat_scRNA-Seq_UMAPs.png"),
  plot = umaps,
  width = 12,
  height = 12,
  dpi = 300,
)

ggsave(
  glue("{workdir}/plots/Seurat_scRNA-Seq_UMAPs.pdf"),
  plot = umaps,
  width = 12,
  height = 12,
  device = "pdf"
)