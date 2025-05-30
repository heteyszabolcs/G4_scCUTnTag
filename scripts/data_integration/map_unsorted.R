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
  "glmGamPoi",
  "ggpubr",
  "RColorBrewer"
)

set.seed(5)

# result folder
result_folder = "../../results/Seurat/unsorted_mousebrain/res0.1/integration/plots/"

# Seurat objects
rna = "../../data/scRNA-Seq/scRNA_Seq-Zeisel_et_al-neuron.Rds"
bartosovic_rna = "../../data/scRNA-Seq/scRNA_Seq-mouse_brain.Rds"
unsorted = readRDS("../../results/Seurat/unsorted_mousebrain/res0.1/outputs/Seurat_object.Rds")
sorted = readRDS("../../results/Seurat/GFP_sorted_mousebrain/res0.8/outputs/Seurat_object.Rds")

# add labels to sorted cluster
sorted@meta.data = sorted@meta.data %>% 
  mutate(seurat_clusters = case_when(seurat_clusters == "0" ~ "OEC/OPC",
                                     seurat_clusters == "1" ~ "OPC",
                                     seurat_clusters == "2" ~ "AST",
                                     seurat_clusters == "3" ~ "OEC"))


rna = readRDS(rna)
all.genes = rownames(rna)
rna = NormalizeData(rna,
                    normalization.method = "LogNormalize",
                    scale.factor = 10000)
rna = FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
rna = ScaleData(rna, features = all.genes)
rna = RunPCA(rna, features = VariableFeatures(object = rna))
rna = RunUMAP(rna, dims = 1:20, return.model = TRUE)

bartosovic_rna = readRDS(bartosovic_rna)
all.genes = rownames(bartosovic_rna)
bartosovic_rna = NormalizeData(bartosovic_rna,
                    normalization.method = "LogNormalize",
                    scale.factor = 10000)
bartosovic_rna = FindVariableFeatures(bartosovic_rna, selection.method = "vst", nfeatures = 2000)
bartosovic_rna = ScaleData(bartosovic_rna, features = all.genes)
bartosovic_rna = RunPCA(bartosovic_rna, features = VariableFeatures(object = bartosovic_rna))
bartosovic_rna = RunUMAP(bartosovic_rna, dims = 1:20, return.model = TRUE)

# integration of unsorted data with neuron scRNA-Seq
# set assays
DefaultAssay(unsorted) = "GA"
DefaultAssay(rna) = "RNA"

common.genes = intersect(rownames(rna), rownames(unsorted))
anchors = FindTransferAnchors(
  reference = rna,
  query = unsorted,
  reduction = 'cca',
  query.assay = 'GA',
  reference.assay = 'RNA',
  k.filter = NA,
  features = common.genes
)

unsorted_neuron = MapQuery(
  anchorset = anchors,
  query = unsorted,
  reference = rna,
  refdata = list(cell_type = "cell_type"),
  reference.reduction = "umap", 
  reduction.model = "umap"
)

p1 = DimPlot(unsorted_neuron, reduction = "ref.umap", group.by = "predicted.cell_type", label = TRUE, label.size = 3, repel = TRUE)
p2 = DimPlot(unsorted_neuron, reduction = "ref.umap", group.by = "seurat_clusters", label = TRUE, label.size = 3, repel = TRUE)
p3 = DimPlot(rna, reduction = "umap", group.by = "cell_type", label = TRUE, label.size = 3, repel = TRUE)
p4 = FeaturePlot(unsorted_neuron, reduction = "ref.umap", features = "predicted.cell_type.score", label = FALSE, label.size = 3, repel = TRUE)
p4

ps = p1 + p2 + p3 + p4
ps

meta = unsorted_neuron@meta.data
high_pred = meta %>% dplyr::filter(predicted.cell_type.score > 0.75)

ggplot(high_pred, aes(y=predicted.cell_type.score, x=seurat_clusters, fill=seurat_clusters)) + 
  geom_bar(stat = "identity")

ggplot(meta, aes(x=predicted.cell_type, y=predicted.cell_type.score, fill=seurat_clusters)) + 
  geom_boxplot() +
  facet_wrap(~predicted.cell_type, scale="free")

ggplot(meta, aes(x=seurat_clusters, y=predicted.cell_type.score, fill=seurat_clusters)) + 
  geom_boxplot() 

# map only unsorted cluster 1
unsorted_cluster1 = subset(x = unsorted_neuron, subset = seurat_clusters == "1")
saveRDS(unsorted_cluster1, "../../results/Seurat/unsorted_mousebrain/res0.1/outputs/Seurat_object_cl1.Rds")

# set assays
DefaultAssay(unsorted_cluster1) = "GA"
DefaultAssay(rna) = "RNA"

common.genes = intersect(rownames(rna), rownames(unsorted_cluster1))
anchors = FindTransferAnchors(
  reference = rna,
  query = unsorted_cluster1,
  reduction = 'cca',
  query.assay = 'GA',
  reference.assay = 'RNA',
  k.filter = NA,
  features = common.genes
)

unsorted_cluster1 <- MapQuery(
  anchorset = anchors,
  query = unsorted_cluster1,
  reference = rna,
  refdata = list(cell_type = "cell_type"),
  reference.reduction = "umap", 
  reduction.model = "umap"
)

p1_cl1 = DimPlot(unsorted_cluster1, reduction = "ref.umap", group.by = "predicted.cell_type", label = TRUE, label.size = 3, repel = TRUE)
p2_cl1 = DimPlot(unsorted_cluster1, reduction = "ref.umap", group.by = "seurat_clusters", label = TRUE, label.size = 3, repel = TRUE)
p3_cl1 = DimPlot(rna, reduction = "umap", group.by = "cell_type", label = TRUE, label.size = 3, repel = TRUE)
p4_cl1 = FeaturePlot(unsorted_cluster1, reduction = "ref.umap", features = "predicted.cell_type.score", label = FALSE, label.size = 3, repel = TRUE)
p4_cl1

ps_cl1 = p1_cl1 + p2_cl1 + p3_cl1 + p4_cl1
ps_cl1


meta = unsorted_cluster1@meta.data
high_pred = meta %>% dplyr::filter(predicted.cell_type.score > 0.5)

ggplot(high_pred, aes(y=predicted.cell_type.score, x=predicted.cell_type, fill=seurat_clusters)) + 
  geom_bar(stat = "identity") 

ggplot(meta, aes(x=predicted.cell_type, y=predicted.cell_type.score, fill=seurat_clusters)) + 
  geom_boxplot() 

# integration of unsorted data with Bartosovic scRNA-Seq
# set assays
DefaultAssay(unsorted) = "GA"
DefaultAssay(bartosovic_rna) = "RNA"

common.genes = intersect(rownames(bartosovic_rna), rownames(unsorted))
anchors = FindTransferAnchors(
  reference = bartosovic_rna,
  query = unsorted,
  reduction = 'cca',
  query.assay = 'GA',
  reference.assay = 'RNA',
  k.filter = NA,
  features = common.genes
)

unsorted_bartosovic = MapQuery(
  anchorset = anchors,
  query = unsorted,
  reference = bartosovic_rna,
  refdata = list(cell_type = "cell_type"),
  reference.reduction = "umap", 
  reduction.model = "umap"
)

# map only unsorted cluster 1
unsorted_cluster0 = subset(x = unsorted_bartosovic, subset = seurat_clusters == "0")
saveRDS(unsorted_cluster0, "../../results/Seurat/unsorted_mousebrain/res0.1/outputs/Seurat_object_cl0.Rds")

# mapping sorted G4 scCUT&Tag to unsorted G4 scCUT&Tag
unsorted = RunTFIDF(unsorted)
unsorted = FindTopFeatures(unsorted, min.cutoff = 'q0')
unsorted = RunSVD(unsorted)
unsorted = RunUMAP(
  unsorted,
  dims = 1:10,
  reduction = 'lsi',
  verbose = FALSE,
  return.model = TRUE
)

# set assays
unsorted[['GA_unsorted']] = unsorted[['GA']]
DefaultAssay(unsorted) = "GA_unsorted"
unsorted[['GA']] = NULL

sorted[['GA_sorted']] = sorted[['GA']]
DefaultAssay(sorted) = "GA_sorted"
sorted[['GA']] = NULL

common.genes = intersect(rownames(sorted), rownames(unsorted))
anchors = FindTransferAnchors(
  reference = unsorted,
  query = sorted,
  reduction = 'cca',
  query.assay = 'GA_sorted',
  reference.assay = 'GA_unsorted',
  k.filter = NA,
  features = common.genes
)

sorted = MapQuery(
  anchorset = anchors,
  query = sorted,
  reference = unsorted,
  refdata = list(seurat_clusters = "seurat_clusters"),
  reference.reduction = "umap", 
  reduction.model = "umap"
)

p1_2 = DimPlot(
  sorted,
  reduction = "umap",
  group.by = "predicted.seurat_clusters",
  label = TRUE,
  label.size = 3,
  repel = TRUE
) +
  labs(
    title = "predicted unsorted labels on GFP+ UMAP",
    x = "UMAP_1",
    y = "UMAP_2",
    fill = NULL
  ) +
  xlim(-5, 8) +
  ylim(-5, 5) +
  scale_color_manual(values = c("#addd8e", "#bcbcbc"))
p1_2

ggsave(
  glue("{result_folder}MapQuery-sorted_to_unsorted-UMAP.pdf"),
  plot = last_plot(),
  device = "pdf",
  width = 7,
  height = 6
)

FeaturePlot(
  sorted,
  reduction = "umap",
  features = "predicted.seurat_clusters.score",
  label = FALSE,
  label.size = 3,
  repel = TRUE
) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  labs(
    title = "predicted unsorted labels on GFP+ UMAP",
    x = "UMAP_1",
    y = "UMAP_2",
    fill = NULL
  ) +
  xlim(-5, 8) +
  ylim(-5, 5) 

ggsave(
  glue("{result_folder}MapQuery-sorted_to_unsorted-predscore_UMAP.pdf"),
  plot = last_plot(),
  device = "pdf",
  width = 7,
  height = 6
)

set3 = brewer.pal(4, "Set3")
cols = c("OEC/OPC"=set3[4], "OPC"=set3[3], "AST"=set3[2], "OEC"=set3[1])
p2_2 = DimPlot(
  sorted,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  label.size = 3,
  repel = TRUE,
  cols = cols
) +
  labs(
    title = "GFP+ UMAP",
    x = "UMAP_1",
    y = "UMAP_2",
    fill = NULL
  ) +
  xlim(-5, 8) +
  ylim(-5, 5)
p2_2

ggsave(
  glue("{result_folder}MapQuery-sorted_to_unsorted-GFPpos_UMAP.pdf"),
  plot = last_plot(),
  device = "pdf",
  width = 7,
  height = 6
)

p3_2 = DimPlot(
  unsorted,
  reduction = "ref.umap",
  group.by = "seurat_clusters",
  label = TRUE,
  label.size = 3,
  repel = TRUE
)
p4_2 = FeaturePlot(
  sorted,
  reduction = "ref.umap",
  features = "predicted.seurat_clusters.score",
  label = FALSE,
  label.size = 3,
  repel = TRUE
)
p4_2

ps_2 = p1_2 + p2_2 + p3_2 + p4_2
ps_2

# visualization of predicted labels
meta = sorted@meta.data
meta %>% group_by(predicted.seurat_clusters) %>% count() %>%
  ggplot(
    .,
    aes(x = predicted.seurat_clusters, y = n, fill =
          predicted.seurat_clusters)
  ) +
  geom_bar(stat="identity", color = "black", show.legend = FALSE) +
  scale_fill_manual(values = c("#addd8e", "#bcbcbc")) +
  ylim(0, 2300) +
  labs(
    title = "",
    x = "predicted cluster",
    y = "number of cells",
    fill = NULL
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 10),
    plot.title = element_text(size = 10),
    axis.title.y = element_text(size = 20, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black")
  ) 

ggsave(
  glue("{result_folder}MapQuery-sorted_to_unsorted-label_quant.pdf"),
  plot = last_plot(),
  device = "pdf",
  width = 7,
  height = 6
)

ggplot(meta,
       aes(x = seurat_clusters, 
           y = predicted.seurat_clusters.score, 
           fill = predicted.seurat_clusters)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#addd8e", "#bcbcbc")) +
  ylim(0.5, 1) +
  labs(
    title = "",
    x = "Seurat cluster (GFP+)",
    y = "prediction score",
    fill = "predicted \nunsorted cluster"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 10),
    plot.title = element_text(size = 10),
    axis.title.y = element_text(size = 20, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.text.x = element_text(size = 20, color = "black", angle = 45, vjust = 0.75, hjust = 0.5),
    axis.text.y = element_text(size = 20, color = "black")
  ) + 
  stat_compare_means(label.y = 1, label = "p.signif")

ggsave(
  glue("{result_folder}MapQuery-sorted_to_unsorted-predscore_distr.pdf"),
  plot = last_plot(),
  device = "pdf",
  width = 7,
  height = 6
)

meta %>% dplyr::filter(predicted.seurat_clusters == 0) %>% pull(predicted.seurat_clusters.score) %>% mean %>% round(2)
meta %>% dplyr::filter(predicted.seurat_clusters == 1) %>% pull(predicted.seurat_clusters.score) %>% mean %>% round(2)


# merging of sorted and unsorted G4 scCUT&Tag datasets
# set assays
# unsorted = readRDS("../../results/Seurat/unsorted_mousebrain/res0.1/outputs/Seurat_object.Rds")
# sorted = readRDS("../../results/Seurat/GFP_sorted_mousebrain/res0.8/outputs/Seurat_object.Rds")
# 
# DefaultAssay(sorted) = "GA"
# DefaultAssay(unsorted) = "GA"
# 
# sorted@meta.data$data_type = "sorted G4 scCnT"
# sorted@meta.data = sorted@meta.data %>% 
#   mutate(custom_col = "sorted clusters")
# 
# 
# unsorted@meta.data$data_type = "unsorted G4 scCnT"
# unsorted@meta.data = unsorted@meta.data %>% mutate(
#   seurat_clusters =
#     case_when(
#       seurat_clusters == "0" ~ "cluster 0 (unsorted)",
#       seurat_clusters == "1" ~ "cluster 1 (unsorted)"
#     )
# ) %>% 
#   mutate(custom_col = seurat_clusters)
# 
# unsorted = FindVariableFeatures(unsorted, selection.method = "vst", nfeatures = 2000)
# genes.use = VariableFeatures(unsorted)
# 
# coembed = merge(x = sorted, y = unsorted)
# coembed = ScaleData(coembed, features = genes.use, do.scale = FALSE)
# coembed = RunPCA(coembed, features = genes.use, verbose = FALSE)
# coembed = RunUMAP(coembed, dims = 1:5)
# 
# cols = c("cluster 0 (unsorted)" = "#fc9272", "cluster 1 (unsorted)" = "#9ecae1",
#          "sorted clusters" = "#31a354")
# coembed_cells = DimPlot(
#   coembed,
#   reduction = "umap",
#   group.by = "custom_col",
#   pt.size = 1,
#   label = FALSE,
#   label.size = 3,
#   repel = TRUE,
#   cols = cols
# ) +
#   labs(
#     title = "sorted + unsorted merged UMAP",
#     x = "UMAP_1",
#     y = "UMAP_2",
#     fill = NULL,
#     cols = cols
#   ) 
# coembed_cells
# 
# ggsave(
#   glue("{result_folder}merged-sorted_unsorted.pdf"),
#   plot = coembed_cells,
#   device = "pdf",
#   width = 7,
#   height = 6
# )
