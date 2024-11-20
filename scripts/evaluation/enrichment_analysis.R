# packages
print("Load R packages")
if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(
  "tidyverse",
  "data.table",
  "ggplot2",
  "glue",
  "topGO",
  "EnsDb.Mmusculus.v79",
  "ComplexHeatmap",
  "circlize",
  "enrichR"
)

#setwd('/crex/proj/snic2020-6-3/SZABOLCS/scG4_study/git_repo/scripts/evaluation')

# source annotation script for mm10
source("../../utils/annotation.R")
source("../../utils/topgo.R")

# export folder
result_folder = "../../results/Seurat/unsorted_mousebrain/"

# unsorted brain scCutnTag Seurat clusters (res 0.1)
peaks_1 = fread("../../results/Seurat/unsorted_mousebrain/res0.1/cluster_spec_peaks/1_peaks.narrowPeak")
peaks_0 = fread("../../results/Seurat/unsorted_mousebrain/res0.1/cluster_spec_peaks/0_peaks.narrowPeak")

cl0 = create_go_matrix(
  create_input(peaks_0),
  colname = "cluster 0"
)

cl1 = create_go_matrix(
  create_input(peaks_1),
  colname = "cluster 1"
)

go_outputs = list(cl0, cl1)

full <- Reduce(function(x, y, ...)
  merge(x, y, all = TRUE, ...),
  go_outputs)

full[is.na(full)] = 1
full = full %>% distinct(Term, .keep_all = TRUE) %>% 
  mutate(`cluster 0` = str_replace(`cluster 0`, "< ", ""), `cluster 1` = str_replace(`cluster 1`, "< ", ""))
rownames(full) = full$Term
full = full %>% dplyr::select(-Term) %>% mutate_if(is.character, as.numeric)
full = as.matrix(full)
full = full[which(rownames(full) != "biological_process"),]

pdf(
  file = glue("{result_folder}topGO_clusterwise_res0.1.pdf"),
  width = 5,
  height = 5
)
print(generate_heatmap(full))
dev.off()

png(
  file = glue("{result_folder}topGO_clusterwise_res0.1.png"),
  width = 11,
  height = 11,
  unit = "cm",
  res = 500
)
print(generate_heatmap(full))
dev.off()

# on FindMarkers outputs
cl1 = fread("../../results/Seurat/unsorted_brain_res0.1-G4markers_logreg-cluster1.tsv")
cl0 = fread("../../results/Seurat/unsorted_brain_res0.1-G4markers_logreg-cluster0.tsv")

cl1 = cl1 %>% dplyr::filter(avg_log2FC > 3)
cl0 = cl0 %>% dplyr::filter(avg_log2FC > 3)

cl1 = create_go_matrix(
  create_input_fam(fc_table = cl1, background = peaks_1),
  colname = "cluster 1"
)
cl0 = create_go_matrix(
  create_input_fam(fc_table = cl0, background = peaks_0),
  colname = "cluster 0"
)

go_outputs = list(cl0, cl1)

full <- Reduce(function(x, y, ...)
  merge(x, y, all = TRUE, ...),
  go_outputs)

full[is.na(full)] = 1
rownames(full) = full$Term
full = full %>% dplyr::select(-Term) %>% mutate_if(is.character, as.numeric)
full = as.matrix(full)
full = full[which(rownames(full) != "biological_process"),]

output = as_tibble(full)
output = output %>% mutate(terms = rownames(full))
write_tsv(output, glue("{result_folder}topGO_input_findmarkers_res0.1.tsv"))

pdf(
  file = glue("{result_folder}topGO_findmarkers_res0.1.pdf"),
  width = 5,
  height = 5
)
print(generate_heatmap_fm(full))
dev.off()

png(
  file = glue("{result_folder}topGO_findmarkers_res0.1.png"),
  width = 11,
  height = 11,
  unit = "cm",
  res = 500
)
print(generate_heatmap_fm(full))
dev.off()

