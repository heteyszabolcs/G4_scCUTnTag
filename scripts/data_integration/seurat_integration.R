print("Load R packages")
if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(
  "Seurat",
  "Signac",
  "glue",
  "ggplot2",
  "dplyr",
  "argparse",
  "tidyverse",
  "ensembldb",
  "EnsDb.Mmusculus.v79",
  "GenomicRanges",
  "cowplot",
  "matrixStats",
  "ggpubr",
  "gridExtra",
  "circlize",
  "ComplexHeatmap"
)

set.seed(5)

# create parser object
parser = ArgumentParser()

parser$add_argument("-w", "--workdir", type = "character", help = "path of working dir")
parser$add_argument("-s", "--seurat_object", type = "character", help = "path to processed scG4 Seurat object with clusters")
parser$add_argument("-r", "--reference", type = "character", help = "select reference scRNA-Seq Seurat object")

args = parser$parse_args()

# add working directory
# work dir for sorted data: "../../results/Seurat/GFP_sorted_mousebrain/res0.8/"
if (all(sapply(args, is.null))) {
  print("Script is running without sbatch.")
  g4 = "../../results/Seurat/GFP_sorted_mousebrain/res0.8/outputs/Seurat_object.Rds"
  workdir = "../../results/Seurat/GFP_sorted_mousebrain/res0.8/integration/"
  rna = "../../data/scRNA-Seq/scRNA_Seq-mouse_brain.Rds"
} else {
  print("Script is running on cluster via bash script.")
  g4 = args$seurat_object
  workdir = args$workdir
  rna = args$reference
}

# make folders for outputs
system(paste0("mkdir -p ", workdir, "/plots"))
system(paste0("mkdir -p ", workdir, "/outputs"))

# G4 - scRNA data integration
# Seurat workflow
print(paste0("Seurat workflow on scRNA-Seq data"))

g4 = readRDS(g4)
rna = readRDS(rna)
all.genes = rownames(rna)
rna = NormalizeData(rna,
                    normalization.method = "LogNormalize",
                    scale.factor = 10000)
rna = FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
rna = ScaleData(rna, features = all.genes)
rna = RunPCA(rna, features = VariableFeatures(object = rna))
rna = RunUMAP(rna, dims = 1:20)

# find all markers in scRNA-Seq data
rna@active.ident = factor(rna$cell_type, levels = unique(rna$cell_type))

markers = FindAllMarkers(rna, min.pct = 0.5)
write_tsv(markers,
          glue("{workdir}/outputs/scRNA-Seq-FindAllMarkers_output.tsv"))
markers.pos = markers[markers$p_val < 0.05 &
                        markers$avg_log2FC > 0.5, ]

# export scRNA Seurat object
saveRDS(rna, glue("{workdir}/outputs/scRNA_Seq_Seurat_object.Rds"))

# set assays
DefaultAssay(g4) = "GA"
DefaultAssay(rna) = "RNA"

common.genes = intersect(rownames(rna), rownames(g4))

# anchor identification between G4 scCnT and scRNA-Seq data sets
transfer.anchors = FindTransferAnchors(
  reference = rna,
  query = g4,
  reduction = 'cca',
  query.assay = 'GA',
  reference.assay = 'RNA',
  k.filter = NA,
  features = common.genes
)

# extract anchor matrix of AnchorSet object
anchor_matrix = as_tibble(transfer.anchors@anchors)
anchor_matrix = anchor_matrix %>%
  mutate(anchor1_barcode = colnames(rna@assays$RNA@counts)[anchor_matrix$cell1]) %>%
  mutate(anchor2_barcode = colnames(g4@assays$GA@counts)[anchor_matrix$cell2])
write_tsv(anchor_matrix, glue("{workdir}/outputs/anchor_matrix.tsv"))

genes.use = VariableFeatures(rna)
refdata = GetAssayData(rna, assay = "RNA", slot = "data")[genes.use, ]

## Transfer data (Seurat)
# imputation - data transfer
imputation = TransferData(
  anchorset = transfer.anchors,
  refdata = refdata,
  weight.reduction = g4[["lsi"]],
  dims = 1:50,
  prediction.assay = TRUE
)

# predict cell labels and add to metadata table
cell_class_pred = TransferData(
  anchorset = transfer.anchors,
  refdata = rna@meta.data$cell_type,
  weight.reduction = g4[["lsi"]],
  dims = 1:50,
  prediction.assay = TRUE
)

predicted_labels = cell_class_pred@data
cell_types = rownames(predicted_labels)
predicted_labels = as.tibble(predicted_labels) %>% mutate(pred_cell_type = cell_types)
predicted_labels = predicted_labels %>% pivot_longer(
  .,
  cols = "AAACGAAAGAAGCCGT-1":"TTTGTGTTCTCGCGTT-1",
  names_to = "cell_id",
  values_to = "pred_max_score"
)
predicted_labels = predicted_labels %>% group_by(cell_id) %>% dplyr::slice(which.max(pred_max_score))

g4_meta = g4@meta.data
g4_meta = g4_meta %>% tibble::rownames_to_column(., var = "barcode") %>%
  inner_join(., predicted_labels, by = c("barcode" = "cell_id"))
barcodes = g4_meta %>% pull(barcode)
g4_meta = g4_meta %>% dplyr::select(-barcode)
rownames(g4_meta) = barcodes
g4@meta.data = g4_meta

# export G4 object with predicted labels
saveRDS(cell_class_pred,
        glue("{workdir}/outputs/g4_cell_label_preds.Rds"))

# export integrated Seurat object
g4[['RNA']] = imputation
saveRDS(g4, glue("{workdir}/outputs/G4_scRNA_integration.Rds"))
