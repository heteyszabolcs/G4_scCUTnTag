print("Load R packages")
if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(
  "data.table",
  "Seurat",
  "rtracklayer",
  "glue",
  "tidyverse",
  "GenomicRanges",
  "cotools",
  "TxDb.Mmusculus.UCSC.mm10.knownGene",
  "ggpubr",
  "ComplexHeatmap",
  "circlize",
  "GenomicFeatures",
  "cotools",
  "bedscout"
)

# helper function for upset plot
source("../../utils/upset_plot.R")

# export folder
result_folder = "../../results/peak_overlaps/"

## mESC-MEF scCut&Tag ###
mesc_mef_cluster0 = fread("../../results/Seurat/mESC_MEF/cluster_spec_peaks/0_peaks.narrowPeak")
mesc_mef_cluster0$type = "cluster0"
mesc_mef_cluster0 = GRanges(
  seqnames = mesc_mef_cluster0$V1,
  ranges = IRanges(
    start = mesc_mef_cluster0$V2,
    end = mesc_mef_cluster0$V3,
    names = mesc_mef_cluster0$type,
  )
)

mesc_mef_cluster1 = fread("../../results/Seurat/mESC_MEF/cluster_spec_peaks/1_peaks.narrowPeak")
mesc_mef_cluster1$type = "cluster1"
mesc_mef_cluster1 = GRanges(
  seqnames = mesc_mef_cluster1$V1,
  ranges = IRanges(
    start = mesc_mef_cluster1$V2,
    end = mesc_mef_cluster1$V3,
    names = mesc_mef_cluster1$type,
  )
)

mef = fread("../../data/bulk_CUTnTag/bulk_G4_CnT_MEF_rep1_R1.mLb.clN_peaks.broadPeak")
mef$type = "MEF"
mef = GRanges(
  seqnames = mef$V1,
  ranges = IRanges(
    start = mef$V2,
    end = mef$V3,
    names = mef$type,
  )
)

mesc = fread("../../data/bulk_CUTnTag/bulk_G4_CnT_mESC_rep1_R1.mLb.clN_peaks.broadPeak")
mesc$type = "mESC"
mesc = GRanges(
  seqnames = mesc$V1,
  ranges = IRanges(
    start = mesc$V2,
    end = mesc$V3,
    names = mesc$type,
  )
)

mef_mesc_upset_input = upset_input(peak_set_1 = mesc_mef_cluster0, 
                               peak_set_2 = mesc_mef_cluster1, 
                               bulk_set_1 = mesc,
                               bulk_set_2 = mef)
mef_mesc_upset_plot = plot_freq_intersect(
  mef_mesc_upset_input,
  .by = "group",
  .levels = c("cluster 0", "cluster 1", "bulk mESC", "bulk MEF"),
  .split = "category",
  .color = "#636363",
  top_n = 10
)
mef_mesc_upset_plot

ggsave(
  glue("{result_folder}upset_plot-mESC_MEF.pdf"),
  plot = mef_mesc_upset_plot,
  width = 6,
  height = 6,
  dpi = 500,
)

## unsorted brain scCut&Tag ###
# unsorted_cluster0 = fread("../../results/Seurat/unsorted_mousebrain/res0.1/cluster_spec_peaks/0_peaks.narrowPeak")
# unsorted_cluster0$type = "cluster0"
# unsorted_cluster0 = GRanges(
#   seqnames = unsorted_cluster0$V1,
#   ranges = IRanges(
#     start = unsorted_cluster0$V2,
#     end = unsorted_cluster0$V3,
#     names = unsorted_cluster0$type,
#   )
# )
# 
# unsorted_cluster1 = fread("../../results/Seurat/unsorted_mousebrain/res0.1/cluster_spec_peaks/1_peaks.narrowPeak")
# unsorted_cluster1$type = "cluster1"
# unsorted_cluster1 = GRanges(
#   seqnames = unsorted_cluster1$V1,
#   ranges = IRanges(
#     start = unsorted_cluster1$V2,
#     end = unsorted_cluster1$V3,
#     names = unsorted_cluster1$type,
#   )
# )
# 
# npc = fread("../../data/bulk_CUTnTag/bulk_CnT_G4_NPC_mm10.bed")
# npc$type = "NPC"
# npc = GRanges(
#   seqnames = npc$V1,
#   ranges = IRanges(
#     start = npc$V2,
#     end = npc$V3,
#     names = npc$type,
#   )
# )
# 
# neuron = fread("../../data/bulk_CUTnTag/bulk_CnT_G4_neuron_mm10.bed")
# neuron$type = "neuron"
# neuron = GRanges(
#   seqnames = neuron$V1,
#   ranges = IRanges(
#     start = neuron$V2,
#     end = neuron$V3,
#     names = neuron$type,
#   )
# )
# 
# npc_upset_input = upset_input(peak_set_1 = unsorted_cluster0, peak_set_2 = unsorted_cluster1, bulk_set = npc)
# npc_upset_plot = plot_freq_intersect(
#   npc_upset_input,
#   .by = "group",
#   .levels = c("cluster 0", "cluster 1", "bulk NPC"),
#   .split = "category",
#   .color = npc_col,
#   top_n = 10
# )
# npc_upset_plot
# 
# neuron_upset_input = upset_input(peak_set_1 = unsorted_cluster0, peak_set_2 = unsorted_cluster1, bulk_set = neuron)
# neuron_upset_plot = plot_freq_intersect(
#   neuron_upset_input,
#   .by = "group",
#   .levels = c("cluster 0", "cluster 1", "bulk neuron"),
#   .split = "category",
#   .color = neuron_col,
#   top_n = 10
# )
# neuron_upset_plot
# 
# ggsave(
#   glue("{result_folder}upset_plot-unsorted_brain-bulkNPC.pdf"),
#   plot = npc_upset_plot,
#   width = 6,
#   height = 6,
#   dpi = 500,
# )
# ggsave(
#   glue("{result_folder}upset_plot-unsorted_brain-bulkNeuron.pdf"),
#   plot = neuron_upset_plot,
#   width = 6,
#   height = 6,
#   device = "pdf"
# )

## Venn diagrams
# mESC-MEF
mef_venn = plot_euler(
  list(mef, mesc_mef_cluster0, mesc_mef_cluster1),
  ignore.strand = TRUE,
  fills = c("#f0f0f0", "#9ecae1", "#fc9272"),
  names = c("MEF", "cluster 0", "cluster 1")
)
ggsave(
  glue("{result_folder}Venn-mESC_MEF-bulkMEF.pdf"),
  plot = mef_venn,
  width = 6,
  height = 6,
  device = "pdf"
)

mesc_venn = plot_euler(
  list(mesc, mesc_mef_cluster0, mesc_mef_cluster1),
  ignore.strand = TRUE,
  fills = c("#f0f0f0", "#9ecae1", "#fc9272"),
  names = c("mESC", "cluster 0", "cluster 1")
)
ggsave(
  glue("{result_folder}Venn-mESC_MEF-bulkmESC.pdf"),
  plot = mesc_venn,
  width = 6,
  height = 6,
  device = "pdf"
)


### jaccard indexes
# mESC-MEF 
peaks = list(mesc_mef_cluster0, mesc_mef_cluster1, mesc, mef)
mesc_mef_jaccards = matrix(NA_real_, length(peaks), length(peaks))
colnames(mesc_mef_jaccards) = c("cluster 0", "cluster 1", "bulk mESC CnT", "bulk MEF CnT")
rownames(mesc_mef_jaccards) = c("cluster 0", "cluster 1", "bulk mESC CnT", "bulk MEF CnT")
for (i in seq(1, ncol(mesc_mef_jaccards))) {
  for (j in seq(1, nrow(mesc_mef_jaccards))) {
    jaccard = genomicCorr.jaccard(peaks[[i]], peaks[[j]])
    mesc_mef_jaccards[i, j] = jaccard
  }
}
is.matrix(mesc_mef_jaccards)

pdf(file = glue("{result_folder}mESC_MEF_res0.1_Jaccard_hm.pdf"),
    width = 6,
    height = 6)
col_fun = colorRamp2(c(0, 0.25, 0.5), c("#9ecae1", "white", "#fc9272"))
Heatmap(
  mesc_mef_jaccards,
  column_title = "",
  row_title = "",
  name = "Jaccard index",
  # row_km = 2,
  # column_km = 1,
  clustering_method_rows = "complete",
  col = col_fun,
  rect_gp = gpar(col = "black", lwd = 0.1),
  #top_annotation = ha,
  show_column_dend = TRUE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_dend = TRUE,
  heatmap_width = unit(8, "cm"),
  heatmap_height = unit(8, "cm"),
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 90,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", mesc_mef_jaccards[i, j]), x, y, gp = gpar(fontsize = 10))
  }
)
dev.off()

# unsorted mouse brain
# peaks = list(unsorted_cluster0,
#              unsorted_cluster1,
#              neuron,
#              npc)
# unsorted_jaccards = matrix(NA_real_, length(peaks), length(peaks))
# colnames(unsorted_jaccards) = c("cluster 0", "cluster 1", "bulk neuron CnT", "bulk NPC CnT")
# rownames(unsorted_jaccards) = c("cluster 0", "cluster 1", "bulk neuron CnT", "bulk NPC CnT")
# for (i in seq(1, ncol(unsorted_jaccards))) {
#   for (j in seq(1, nrow(unsorted_jaccards))) {
#     jaccard = genomicCorr.jaccard(peaks[[i]], peaks[[j]])
#     unsorted_jaccards[i, j] = jaccard
#   }
# }
# is.matrix(unsorted_jaccards)
# 
# pdf(
#   file = glue("{result_folder}unsorted_brain_res0.1_Jaccard_hm.pdf"),
#   width = 6,
#   height = 6
# )
# col_fun = colorRamp2(c(0, 0.25, 0.5), c("#9ecae1", "white", "#fc9272"))
# Heatmap(
#   unsorted_jaccards,
#   column_title = "",
#   row_title = "",
#   name = "jaccard index",
#   # row_km = 2,
#   # column_km = 1,
#   clustering_method_rows = "complete",
#   col = col_fun,
#   rect_gp = gpar(col = "black", lwd = 0.1),
#   #top_annotation = ha,
#   show_column_dend = TRUE,
#   cluster_columns = TRUE,
#   cluster_rows = TRUE,
#   show_row_dend = TRUE,
#   heatmap_width = unit(8, "cm"),
#   heatmap_height = unit(8, "cm"),
#   row_names_gp = gpar(fontsize = 10),
#   column_names_gp = gpar(fontsize = 10),
#   column_names_rot = 90,
#   cell_fun = function(j, i, x, y, width, height, fill) {
#     grid.text(sprintf("%.2f", unsorted_jaccards[i, j]), x, y, gp = gpar(fontsize = 10))
#   }
# )
# dev.off()

