if (!require("pacman"))
  install.packages("pacman")
pacman::p_load("Seurat", "glue", "tidyverse", "data.table"
)

# make folders for outputs
system(paste0("mkdir ", "../../results/scBridge/input/"))

# scBridge needs:
# A: gene activity scores of scCut&Tag
# B: normalized scRNA-Seq counts
# keep those genes that overlap between A and B.

rna_path = "../../data/scRNA-Seq/scRNA_Seq-mouse_brain.Rds"
g4_path = "../../results/Seurat/GFP_sorted_mousebrain/res0.8/outputs/Seurat_object.Rds"
result_folder = "../../results/scBridge/"

# scRNA-Seq input - Bartosovic et al.
rna_bartosovic = readRDS(rna_path)
meta = rna_bartosovic@meta.data
rna_bartosovic = rna_bartosovic@assays$RNA@data
rna_bartosovic = as.data.frame(rna_bartosovic)

colnames(rna_bartosovic) = unname(sapply(colnames(rna_bartosovic), function(x) {
  strsplit(x, "_1")[[1]][1]
}))

# G4 input
g4 = readRDS(g4_path)
ga = g4@assays$GA$counts
ga = as.data.frame(ga)

# intersect
int = intersect(rownames(ga), rownames(rna_bartosovic))
rna_bartosovic = rna_bartosovic[int, ]
ga = ga[int, ]

# export
write.csv(
  rna_bartosovic,
  file = glue("{result_folder}GSE163484_Bartosovic_et_al_counts.csv"),
  quote = FALSE
)
cell_ids = unname(sapply(rownames(meta), function(x) {
  strsplit(x, "_1")[[1]][1]
}))
annot = tibble(cell_id = cell_ids, CellType = meta$cell_type)
write_csv(annot,
          file = glue("{result_folder}GSE163484_Bartosovic_et_al_annot.csv"))
write.csv(ga,
          file = glue("{result_folder}GFPsorted-gene_activity_scores.csv"),
          quote = FALSE)
