print("Load R packages")
if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(
  "Seurat",
  "Signac",
  "glue",
  "tidyverse",
  "rtracklayer",
  "EnsDb.Mmusculus.v79",
  "GenomicFeatures",
  "GenomeInfoDb",
  "cowplot",
  "Matrix",
  "ggpubr",
  "DoubletFinder",
  "argparse",
  "hrbrthemes"
)

set.seed(5)

# create parser object
parser = ArgumentParser()

parser$add_argument("-c",
                    "--cell_ranger_output",
                    type = "character",
                    help = "path to output folder of Cell Ranger")
parser$add_argument("-w", "--workdir", type = "character", help = "path of working dir")
parser$add_argument("-r", "--res", type = "numeric", help = "cluster resolution")
parser$add_argument("-s", "--sample", type = "character", help = "sample: mESC_MEF, unsorted, sorted")

args = parser$parse_args()

# add working directory
# e.g. "../../results/Seurat/mESC_MEF/" or "../../results/Seurat/unsorted_mousebrain/res0.8/"
# add path to CellRanger output folder
# e.g. "../../data/CellRanger/mESC_MEF/"
# cluster resolution (0.1 or 0.8 were used)
if (all(sapply(args, is.null))) {
  print("Script is running without sbatch.")
  res = 0.1
  workdir = "../../results/Seurat/mESC_MEF/"
  cell_ranger = "../../data/CellRanger/mESC_MEF/"
  sample = "mESC_MEF"
} else {
  print("Script is running on cluster via bash script.")
  cell_ranger = args$cell_ranger_output
  workdir = args$workdir
  res = args$res
  sample = args$sample
}

print(paste0("Seurat cluster resolution: ", as.character(res)))

# make folders for outputs
system(paste0("mkdir -p ", workdir))
system(paste0("mkdir -p ", workdir, "/plots"))
system(paste0("mkdir -p ", workdir, "/outputs"))

# read into Seurat
counts = Read10X_h5(filename = glue("{cell_ranger}filtered_peak_bc_matrix.h5"))
metadata = read.csv(
  file = glue("{cell_ranger}singlecell.csv"),
  header = TRUE,
  row.names = 1
)

## Seurat workflow
print(paste0("Run Seurat workflow on ", cell_ranger))
chrom_assay = CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = glue("{cell_ranger}fragments.tsv.gz"),
  min.cells = 10,
  min.features = 200
)

seurat = CreateSeuratObject(counts = chrom_assay,
                            assay = "peaks",
                            meta.data = metadata)

# normalization
seurat = RunTFIDF(seurat)
seurat = FindTopFeatures(seurat, min.cutoff = 'q0')
seurat = RunSVD(seurat)
depthcor = DepthCor(seurat) # 1st component highly correlate with sequencing depth

# Non-linear dimension reduction and clustering
seurat = RunUMAP(object = seurat,
                 reduction = 'lsi',
                 dims = 2:30)

# find doublets using DoubletFinder
doublet_test = Seurat::Read10X_h5(
  filename = glue("{cell_ranger}/filtered_peak_bc_matrix.h5"),
  use.names = T
)
doublet_test = CreateSeuratObject(doublet_test, project = "G4 scC&T")
doublet_test = NormalizeData(doublet_test)
doublet_test = FindVariableFeatures(doublet_test, verbose = F)
doublet_test = ScaleData(
  doublet_test,
  vars.to.regress = c("nFeature_RNA", "percent_mito"),
  verbose = F
)
doublet_test = RunPCA(doublet_test, verbose = F, npcs = 20)
doublet_test = RunUMAP(doublet_test, dims = 1:10, verbose = F)

nExp = round(ncol(doublet_test) * 0.04)  # expect 4% doublets
# if doubletFinder function does not work due to Seurat versions:
# remotes::install_github(“egeulgen/DoubletFinder”)
doublet_test <- doubletFinder_v3(
  doublet_test,
  pN = 0.25,
  pK = 0.09,
  nExp = nExp,
  PCs = 1:10
)
singlets = rownames(doublet_test@meta.data[which(doublet_test@meta.data[, 5] == "Singlet"), ])

# keep only singlets!
meta = seurat@meta.data
seurat@meta.data = meta %>% mutate(doublet_test = ifelse(rownames(meta) %in% singlets, "Singlet", NA_character_))
seurat = subset(x = seurat, subset = doublet_test == "Singlet")

# give EnsDB mm10 annotation
annotations = GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) = 'UCSC'
Annotation(seurat) = annotations

# clustering (community detection)
seurat = FindNeighbors(object = seurat,
                       reduction = 'lsi',
                       dims = 2:30)
seurat = FindClusters(
  object = seurat,
  verbose = FALSE,
  resolution = res,
  algorithm = 3
)

# in the case of mouse brain data
# removing clusters with very low cell number
if (res == 0.8) {
  seurat = subset(x = seurat,
                  idents = 4,
                  invert = TRUE)
}

# QC violin plots (without correction)
nF_violin = VlnPlot(seurat,
                    group.by = "seurat_clusters",
                    features = "nFeature_peaks",
                    pt.size = 0.1) +
  #scale_fill_manual(values = c("#addd8e", "#bdbdbd", "#addd8e")) +
  ggtitle("nFeature (peaks)") +
  xlab("cluster") +
  ylim(0, 5000) +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(
      size = 25,
      color = "black",
      angle = 0
    ),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoLegend()
nC_violin = VlnPlot(seurat,
                    group.by = "seurat_clusters",
                    features = "nCount_peaks",
                    pt.size = 0.1) +
  #scale_fill_manual(values = c("#addd8e", "#bdbdbd", "#addd8e")) +
  ggtitle("nCount (peaks)") +
  xlab("cluster") +
  ylim(0, 20000) +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(
      size = 25,
      color = "black",
      angle = 0
    ),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoLegend()
TSS_violin = VlnPlot(seurat,
                     group.by = "seurat_clusters",
                     features = "TSS_fragments",
                     pt.size = 0.1) +
  #scale_fill_manual(values = c("#addd8e", "#bdbdbd", "#addd8e")) +
  ggtitle("TSS fragments") +
  xlab("cluster") +
  ylim(0, 10000) +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(
      size = 25,
      color = "black",
      angle = 0
    ),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoLegend()
mito_violin = VlnPlot(seurat,
                      group.by = "seurat_clusters",
                      features = "mitochondrial",
                      pt.size = 0.1) +
  #scale_fill_manual(values = c("#addd8e", "#bdbdbd", "#addd8e")) +
  ggtitle("mitochondrial fragments") +
  xlab("cluster") +
  ylim(0, 30000) +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(
      size = 25,
      color = "black",
      angle = 0
    ),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoLegend()
qc_violins = ggarrange(nF_violin, nC_violin, TSS_violin, mito_violin)

ggsave(
  glue("{workdir}/plots/QC_violins.pdf"),
  plot = qc_violins,
  width = 12,
  height = 10,
  dpi = 300,
)

ggsave(
  glue("{workdir}/plots/QC_violins.pdf"),
  plot = qc_violins,
  device = "pdf",
  width = 12,
  height = 10,
  dpi = 300,
)

# exclude outlier cells
clusters = levels(seurat@meta.data$seurat_clusters)
if (res == 0.1) {
  if (2 %in% clusters) {
    seurat = subset(x = seurat,
                    idents = 2,
                    invert = TRUE)
  }
}

# marker analysis
# by logistic regression with total number of fragments as a latent variable
print(paste0("Marker analysis"))
markers = FindAllMarkers(seurat, test.use = "LR", latent.vars = "peak_region_fragments")
markers = markers %>% dplyr::filter(p_val_adj < 0.05)

markers = markers %>% separate(
  gene,
  sep = "-",
  into = c("chr", "start", "end"),
  remove = TRUE
)
write_tsv(markers,
          glue("{workdir}/outputs/FindAllMarkers_logreg_output.tsv"))
bed = markers %>% dplyr::select("chr", "start", "end")
write_tsv(bed,
          glue("{workdir}/outputs/FindAllMarkers_logreg_output.bed"),
          col_names = FALSE)

# by wilcoxon
wilcoxon = FindAllMarkers(seurat, test.use = "wilcox")
wilcoxon = wilcoxon %>% dplyr::filter(p_val_adj < 0.05)
write_tsv(wilcoxon,
          glue("{workdir}/outputs/FindAllMarkers_wilcoxon_output.tsv"))

# cell numbers per cluster
clusters = unique(seurat$seurat_clusters)
for (i in clusters) {
  length = seurat@meta.data %>% dplyr::filter(seurat_clusters == as.numeric(i)) %>% rownames %>% length
  print(glue("Number of cells in cluster {i}: {length}"))
}

for (i in clusters) {
  median = seurat@meta.data %>% dplyr::filter(seurat_clusters == i) %>% pull(nFeature_peaks) %>% median
  print(glue("Median read count of features in cluster {i}: {median}"))
}

# plot UMAPs
if (res == 0.1 & sample == "mESC_MEF") {
  cols = c("0" = "#9ecae1", "1" = "#fc9272")
  
  dim_res0.1 = DimPlot(
    object = seurat,
    label = TRUE,
    pt.size = 2,
    label.size = 7,
    repel = TRUE,
    raster = TRUE
  ) +
    xlim(-10, 10) +
    ylim(-10, 10) +
    scale_colour_manual(values = cols, breaks = c("0", "1"), labels = c("MEF", "mESC")) +
                        ggtitle("") +
                          theme(
                            text = element_text(size = 25),
                            plot.title = element_text(size = 20),
                            axis.text.x = element_text(size = 25, color = "black"),
                            axis.text.y = element_text(size = 25, color = "black")
                          )
                        
                        ggsave(
                          glue("{workdir}/plots/Seurat_UMAP_res0.1.png"),
                          plot = dim_res0.1,
                          width = 10,
                          height = 10,
                          dpi = 300,
                        )
                        
                        ggsave(
                          glue("{workdir}/plots/Seurat_UMAP_res0.1.pdf"),
                          plot = dim_res0.1,
                          width = 10,
                          height = 10,
                          dpi = 300,
                        )
} 

if (res == 0.1 & sample == "unsorted") {
  cols = c("0" = "#addd8e", "1" = "#bdbdbd")
  
  dim_res0.1 = DimPlot(
    object = seurat,
    label = TRUE,
    pt.size = 2,
    label.size = 7,
    repel = TRUE,
    raster = TRUE
  ) +
    xlim(-10, 10) +
    ylim(-10, 10) +
    scale_colour_manual(values = cols, breaks = c("0", "1"), labels = c("cluster 0", "cluster 1")) +
    ggtitle("") +
    theme(
      text = element_text(size = 25),
      plot.title = element_text(size = 20),
      axis.text.x = element_text(size = 25, color = "black"),
      axis.text.y = element_text(size = 25, color = "black")
    )
  
  ggsave(
    glue("{workdir}/plots/Seurat_UMAP_res0.1.png"),
    plot = dim_res0.1,
    width = 10,
    height = 10,
    dpi = 300,
  )
  
  ggsave(
    glue("{workdir}/plots/Seurat_UMAP_res0.1.pdf"),
    plot = dim_res0.1,
    width = 10,
    height = 10,
    dpi = 300,
  )
} else {
  dim = DimPlot(
    object = seurat,
    label = TRUE,
    pt.size = 1,
    label.size = 7,
    repel = TRUE,
    order = c(0, 1, 2, 3)
  ) +
    NoLegend() +
    scale_color_brewer(palette = "Set3") +
    xlim(-10, 10) +
    ylim(-10, 10) +
    ggtitle("") +
    theme(
      text = element_text(size = 25),
      plot.title = element_text(size = 20),
      axis.text.x = element_text(size = 25, color = "black"),
      axis.text.y = element_text(size = 25, color = "black")
    )
  
  dim_blank = DimPlot(
    object = seurat,
    label = TRUE,
    pt.size = 1,
    label.size = 7,
    repel = TRUE,
    order = c(0, 1, 2, 3)
  ) +
    NoLegend() +
    NoAxes() +
    scale_color_brewer(palette = "Set3") +
    xlim(-10, 10) +
    ylim(-10, 10)
  
  ggsave(
    glue("{workdir}/plots/Seurat_UMAP.png"),
    plot = dim,
    width = 10,
    height = 10,
    dpi = 300
  )
  
  ggsave(
    glue("{workdir}/plots/Seurat_UMAP.pdf"),
    plot = dim,
    device = "pdf",
    width = 10,
    height = 10
  )
  
  ggsave(
    glue("{workdir}/plots/Seurat_UMAP_blankedplot.pdf"),
    plot = dim_blank,
    width = 10,
    height = 10,
    device = "pdf"
  )
  
  ggsave(
    glue("{workdir}/plots/Seurat_UMAP_blankedplot.png"),
    plot = dim_blank,
    width = 10,
    height = 10,
    dpi = 300
  )
  
}

# indicate doublets
print("UMAP of DoubletFinder")
DF.name = colnames(doublet_test@meta.data)[grepl("DF.classification", colnames(doublet_test@meta.data))]
doublets = DimPlot(doublet_test,
                   group.by = DF.name,
                   pt.size = 0.3,
                   order = "Doublet") +
  scale_color_manual(values = c("grey", "#de2d26")) +
  xlim(-10, 10) +
  ylim(-10, 10) +
  ggtitle("Doublet finder")

doublet_ctrl = DimPlot(doublet_test, group.by = "orig.ident", pt.size = 0.3) +
  scale_color_manual(values = "#de2d26") +
  xlim(-10, 10) +
  ylim(-10, 10) +
  ggtitle("")

ggsave(
  glue("{workdir}/plots/Seurat_doublets.png"),
  plot = doublets,
  width = 5,
  height = 5,
  dpi = 300
)

ggsave(
  glue("{workdir}/plots/Seurat_doublets.pdf"),
  device = "pdf",
  plot = doublets,
  width = 5,
  height = 5
)

ggsave(
  glue("{workdir}/plots/Seurat_doublets_ctrl.png"),
  plot = doublet_ctrl,
  width = 5,
  height = 5,
  dpi = 300
)

ggsave(
  glue("{workdir}/plots/Seurat_doublets_ctrl.pdf"),
  device = "pdf",
  plot = doublet_ctrl,
  width = 5,
  height = 5
)

# Create a gene activity matrix
# GeneActivity scores: counts over 2 kb-upstream region and gene body
gene.activities = GeneActivity(seurat)

seurat[['GA']] = CreateAssayObject(counts = gene.activities)
seurat = NormalizeData(
  object = seurat,
  assay = 'GA',
  normalization.method = 'LogNormalize',
  scale.factor = median(seurat$nCount_peaks)
)

# find marker genes at the level of GA
ga_markers = FindAllMarkers(seurat, only.pos = TRUE, assay = "GA")
write_tsv(ga_markers,
          glue("{workdir}/outputs/FindAllMarkers_GA_output.tsv"))

ga_markers_top = ga_markers %>% group_by(cluster) %>% top_n(wt = avg_log2FC, n = 5)
ga_markers_hm = DoHeatmap(
  object = seurat,
  features = ga_markers_top$gene,
  slot = 'data',
  raster = TRUE,
  assay = "GA"
) +
  scale_fill_gradient2(low = "white",
                       mid = "white",
                       high = "red") +
  labs(title = " ", fill = "gene act. score") +
  theme(plot.title = element_text(size = 8),
        axis.text.y = element_text(size = 13, color = "black"))

ggsave(
  glue("{workdir}/plots/GA_heatmap.pdf"),
  plot = ga_markers_hm,
  device = "pdf",
  width = 10,
  height = 3
)

# export Rds
saveRDS(seurat, glue("{workdir}/outputs/Seurat_object.Rds"))
