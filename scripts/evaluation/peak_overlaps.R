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
unsorted_cluster0 = fread("../../results/Seurat/unsorted_mousebrain/res0.1/cluster_spec_peaks/0_peaks.narrowPeak")
unsorted_cluster0$type = "cluster0"
unsorted_cluster0 = GRanges(
  seqnames = unsorted_cluster0$V1,
  ranges = IRanges(
    start = unsorted_cluster0$V2,
    end = unsorted_cluster0$V3,
    names = unsorted_cluster0$type,
  )
)

unsorted_cluster1 = fread("../../results/Seurat/unsorted_mousebrain/res0.1/cluster_spec_peaks/1_peaks.narrowPeak")
unsorted_cluster1$type = "cluster0"
unsorted_cluster1 = GRanges(
  seqnames = unsorted_cluster1$V1,
  ranges = IRanges(
    start = unsorted_cluster1$V2,
    end = unsorted_cluster1$V3,
    names = unsorted_cluster1$type,
  )
)

## sorted brain scCut&Tag ###
sorted = fread("../../data/CellRanger/GFP_sorted_mousebrain/peaks.bed")
sorted$type = "GFP+ G4 scCUT&Tag"
sorted = GRanges(
  seqnames = sorted$V1,
  ranges = IRanges(
    start = sorted$V2,
    end = sorted$V3,
    names = sorted$type,
  )
)

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

# unsorted vs. sorted
unsorted_sorted_venn = plot_euler(
  list(sorted, unsorted_cluster0, unsorted_cluster1),
  ignore.strand = TRUE,
  fills = c("#ffffff", "#addd8e", "#bcbcbc"),
  names = c("sorted", "unsorted cl. 0", "cl. 1")
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