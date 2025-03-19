if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(
  "EnsDb.Mmusculus.v75",
  "Seurat",
  "Signac",
  "ggpubr",
  "glue",
  "RColorBrewer",
  "tidyverse"
)

# function for creating quality control violin plots
qc = function(filtered_peak_mat = "../../data/CellRanger/unsorted_mousebrain/filtered_peak_bc_matrix.h5",
              metadata = "../../data/CellRanger/unsorted_mousebrain/singlecell.csv",
              fragments = "../../data/CellRanger/unsorted_mousebrain/fragments.tsv.gz",
              ensdb = EnsDb.Mmusculus.v75,
              label = "unsorted",
              res = 0.1) {
  counts <- Read10X_h5(filename = filtered_peak_mat)
  metadata <- read.csv(file = metadata,
                       header = TRUE,
                       row.names = 1)
  chrom_assay = CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    genome = "mm10",
    fragments = fragments,
    min.cells = 10,
    min.features = 200
  )
  
  seurat = CreateSeuratObject(counts = chrom_assay,
                              assay = "peaks",
                              meta.data = metadata)
  
  # nucleosome signal score
  seurat = NucleosomeSignal(object = seurat)
  
  # compute TSS enrichment score per cell
  annotations = GetGRangesFromEnsDb(ensdb = ensdb)
  seqlevels(annotations) = paste0('chr', seqlevels(annotations))
  genome(annotations) = "mm10"
  
  Annotation(seurat) = annotations
  seurat = TSSEnrichment(object = seurat, fast = FALSE)
  
  # add blacklist ratio and fraction of reads in peaks (FRIP)
  seurat$pct_reads_in_peaks = seurat$peak_region_fragments / seurat$passed_filters * 100
  seurat$blacklist_ratio = seurat$blacklist_region_fragments / seurat$peak_region_fragments
  
  DensityScatter(
    seurat,
    x = 'nCount_peaks',
    y = 'TSS.enrichment',
    log_x = TRUE,
    quantiles = TRUE
  )
  ggsave(
    glue(
      "../../results/Seurat/{label}_res{as.character(res)}_density_scatter.pdf"
    ),
    plot = last_plot(),
    width = 6,
    height = 4,
    dpi = 300,
  )
  
  ## processing
  # normalization
  seurat = RunTFIDF(seurat)
  seurat = FindTopFeatures(seurat, min.cutoff = 'q0')
  seurat = RunSVD(seurat)
  
  # Non-linear dimension reduction and clustering
  seurat = RunUMAP(object = seurat,
                   reduction = 'lsi',
                   dims = 2:30)
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
  
  
  seurat$high.tss = ifelse(seurat$TSS.enrichment > 3, 'High', 'Low')
  TSSPlot(seurat, group.by = 'high.tss') + NoLegend()
  
  numb_low_tss = seurat@meta.data %>% dplyr::filter(high.tss == "Low") %>% rownames %>% length
  print(glue("Number of cells with low TSS (TSS threshold: 3): {numb_low_tss}"))
  numb_high_tss = seurat@meta.data %>% dplyr::filter(high.tss == "High") %>% rownames %>% length
  print(glue(
    "Number of cells with high TSS (TSS threshold: 3): {numb_high_tss}"
  ))
  
  print(seurat@meta.data %>% count(seurat_clusters) %>% arrange(desc(n)))
  
  cluster_numb = length(unique(seurat@meta.data$seurat_clusters))
  
  tss_enrich = VlnPlot(
    object = seurat,
    features = 'TSS.enrichment',
    group.by = "seurat_clusters",
    pt.size = 0.1,
    alpha = 0.3,
    ncol = 1,
    cols = if (res != 0.1)
      brewer.pal(cluster_numb, "Set3")
    else
      c("#9ecae1", "#fc9272", "#646a6d")
  ) +
    xlab("cluster") +
    ylab("TSS enrichment") +
    theme(
      text = element_text(size = 15),
      plot.title = element_text(size = 15),
      axis.text.x = element_text(
        size = 25,
        color = "black",
        angle = 0,
        hjust = 0.5
      ),
      axis.text.y = element_text(size = 15, color = "black")
    ) +
    NoLegend()
  nucl_signal = VlnPlot(
    object = seurat,
    features = 'nucleosome_signal',
    group.by = "seurat_clusters",
    pt.size = 0.1,
    alpha = 0.3,
    ncol = 1,
    cols = if (res != 0.1)
      brewer.pal(cluster_numb, "Set3")
    else
      c("#9ecae1", "#fc9272", "#646a6d")
  ) +
    xlab("cluster") +
    ylab("nucleosome signal") +
    theme(
      text = element_text(size = 15),
      plot.title = element_text(size = 15),
      axis.text.x = element_text(
        size = 25,
        color = "black",
        angle = 0,
        hjust = 0.5
      ),
      axis.text.y = element_text(size = 15, color = "black")
    ) +
    NoLegend()
  pct_in_reads = VlnPlot(
    object = seurat,
    features = 'pct_reads_in_peaks',
    group.by = "seurat_clusters",
    pt.size = 0.1,
    alpha = 0.3,
    ncol = 1,
    cols = if (res != 0.1)
      brewer.pal(cluster_numb, "Set3")
    else
      c("#9ecae1", "#fc9272", "#646a6d")
  ) +
    xlab("cluster") +
    ylab("pct. of reads in peaks") +
    theme(
      text = element_text(size = 15),
      plot.title = element_text(size = 15),
      axis.text.x = element_text(
        size = 15,
        color = "black",
        angle = 0,
        hjust = 0.5
      ),
      axis.text.y = element_text(size = 15, color = "black")
    ) +
    NoLegend()
  
  seurat@meta.data = seurat@meta.data %>% mutate(log10_nFeature_peaks = log10(nFeature_peaks))
  nF_violin = VlnPlot(
    seurat,
    group.by = "seurat_clusters",
    features = "log10_nFeature_peaks",
    pt.size = 0.1,
    alpha = 0.3,
    ncol = 1,
    cols = if (res != 0.1)
      brewer.pal(cluster_numb, "Set3")
    else
      c("#9ecae1", "#fc9272", "#646a6d")
  ) +
    xlab("cluster") +
    ylab("number of detected genes [log10]") +
    theme(
      text = element_text(size = 15),
      plot.title = element_text(size = 15),
      axis.text.x = element_text(
        size = 15,
        color = "black",
        angle = 0,
        hjust = 0.5
      ),
      axis.text.y = element_text(size = 15, color = "black")
    ) +
    NoLegend()

  seurat@meta.data = seurat@meta.data %>% mutate(log10_TSS_fragments = log10(TSS_fragments))
  TSS_violin = VlnPlot(
    seurat,
    group.by = "seurat_clusters",
    features = "log10_TSS_fragments",
    pt.size = 0.1,
    alpha = 0.3,
    ncol = 1,
    cols = if (res != 0.1)
      brewer.pal(cluster_numb, "Set3")
    else
      c("#9ecae1", "#fc9272", "#646a6d")
  ) +
    xlab("cluster") +
    ylab("TSS fragment [log10]") +
    theme(
      text = element_text(size = 15),
      plot.title = element_text(size = 15),
      axis.text.x = element_text(
        size = 15,
        color = "black",
        angle = 0,
        hjust = 0.5
      ),
      axis.text.y = element_text(size = 15, color = "black")
    )
  
  seurat@meta.data = seurat@meta.data %>% mutate(log10_mitochondrial = log10(mitochondrial+0.1))
  mito_violin = VlnPlot(
    seurat,
    group.by = "seurat_clusters",
    features = "log10_mitochondrial",
    pt.size = 0.1,
    alpha = 0.3,
    ncol = 1,
    cols = if (res != 0.1)
      brewer.pal(cluster_numb, "Set3")
    else
      c("#9ecae1", "#fc9272", "#646a6d")
  ) +
    xlab("cluster") +
    ylab("mitochondrial fragment [log10]") +
    theme(
      text = element_text(size = 15),
      plot.title = element_text(size = 15),
      axis.text.x = element_text(
        size = 15,
        color = "black",
        angle = 0,
        hjust = 0.5
      ),
      axis.text.y = element_text(size = 15, color = "black")
    )
  
  seurat@meta.data = seurat@meta.data %>% mutate(log10_nCount_peaks = log10(nCount_peaks))
  log10_ncount_violin = VlnPlot(
    seurat,
    group.by = "seurat_clusters",
    features = "log10_nCount_peaks",
    pt.size = 0.1,
    alpha = 0.3,
    ncol = 1,
    cols = if (res != 0.1)
      brewer.pal(cluster_numb, "Set3")
    else
      c("#9ecae1", "#fc9272", "#646a6d")
  ) +
    xlab("cluster") +
    ylab("reads per cell [log10]") +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(
      size = 15,
      color = "black",
      angle = 0,
      hjust = 0.5
    ),
    axis.text.y = element_text(size = 15, color = "black")
  )
  
  qc_violins = ggarrange(
    tss_enrich,
    nucl_signal,
    pct_in_reads,
    nF_violin,
    log10_ncount_violin,
    TSS_violin,
    mito_violin
  )
  
  ggsave(
    glue(
      "../../results/Seurat/{label}_res{as.character(res)}_quality_plots.pdf"
    ),
    plot = qc_violins,
    width = 14,
    height = 12,
    dpi = 300,
  )
  
  return(qc_violins)
  
}

## run
# mouse brain data
qc(
  filtered_peak_mat = "../../data/CellRanger/GFP_sorted_mousebrain/filtered_peak_bc_matrix.h5",
  metadata = "../../data/CellRanger/GFP_sorted_mousebrain/singlecell.csv",
  fragments = "../../data/CellRanger/GFP_sorted_mousebrain/fragments.tsv.gz",
  ensdb = EnsDb.Mmusculus.v75,
  label = "sorted",
  res = 0.8
)

qc(
  filtered_peak_mat = "../../data/CellRanger/unsorted_mousebrain/filtered_peak_bc_matrix.h5",
  metadata = "../../data/CellRanger/unsorted_mousebrain/singlecell.csv",
  fragments = "../../data/CellRanger/unsorted_mousebrain/fragments.tsv.gz",
  ensdb = EnsDb.Mmusculus.v75,
  label = "unsorted",
  res = 0.1
)

qc(
  filtered_peak_mat = "../../data/CellRanger/unsorted_mousebrain/filtered_peak_bc_matrix.h5",
  metadata = "../../data/CellRanger/unsorted_mousebrain/singlecell.csv",
  fragments = "../../data/CellRanger/unsorted_mousebrain/fragments.tsv.gz",
  ensdb = EnsDb.Mmusculus.v75,
  label = "unsorted",
  res = 0.8
)

# mESC-MEF data
# # colors for mESC-MEF plots
# else c("#addd8e", "#fc9272", "#636363")
qc(
  filtered_peak_mat = "../../data/CellRanger/mESC_MEF/filtered_peak_bc_matrix.h5",
  metadata = "../../data/CellRanger/mESC_MEF/singlecell.csv",
  fragments = "../../data/CellRanger/mESC_MEF/fragments.tsv.gz",
  ensdb = EnsDb.Mmusculus.v75,
  label = "mESC-MEF",
  res = 0.1
)
