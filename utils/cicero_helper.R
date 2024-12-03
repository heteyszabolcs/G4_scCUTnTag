suppressMessages({
  library("cicero")
  library("Seurat")
  library("Signac")
  library("ggplot2")
  library("monocle3")
  library("data.table")
  library("tidyverse")
  library("glue")
})

# helper function for cicero analysis
cicero_function = function(seurat_object, label, output_path) {
  input_cds = monocle3::new_cell_data_set(
    expression_data =  seurat_object[['peaks']]@counts,
    cell_metadata = seurat_object@meta.data,
    gene_metadata = data.frame(
      'gene_short_name' = rownames(seurat_object[['peaks']]),
      row.names = rownames(seurat_object[['peaks']])
    )
  )
  input_cds = detect_genes(input_cds)
  input_cds = input_cds[Matrix::rowSums(exprs(input_cds)) != 0, ]
  
  input_cds = detect_genes(input_cds)
  input_cds = estimate_size_factors(input_cds)
  input_cds = preprocess_cds(input_cds, method = "LSI")
  
  # dimension reduction using LSI followed by UMAP
  input_cds = reduce_dimension(input_cds,
                               reduction_method = 'UMAP',
                               preprocess_method = "LSI")
  
  
  umap_coords = reducedDims(input_cds)$UMAP
  cicero_cds = make_cicero_cds(input_cds, reduced_coordinates = umap_coords)
  
  # estimate the co-G4 sites in the genome in order to predict cis-regulatory interactions
  print("Run cicero...")
  conns = run_cicero(cicero_cds, mm10.chromsize)
  save(conns, file = glue("{output_path}cicero_GFPsorted-{label}.Rds"))
  
  # finding cis-Coaccessibility networks (CCANS), coaccess_cutoff_override > 0.1 strict
  CCAN_assigns = generate_ccans(conns, coaccess_cutoff_override = 0.1)
  save(CCAN_assigns, file = glue("{output_path}cicero_GFPsorted_coG4networks-{label}.Rds"))
  
  return(CCAN_assigns)
}

make_coverage_plot = function(region,
                              distance = 5000,
                              seurat_object = g4,
                              result_folder = "../../results/cicero/") {
  chr = strsplit(region, "-")[[1]][1]
  start = as.numeric(strsplit(region, "-")[[1]][2]) - distance
  end = as.numeric(strsplit(region, "-")[[1]][3]) + distance
  p = CoveragePlot(
    seurat_object,
    region = glue("{chr}-{start}-{end}"),
    annotation = TRUE,
    show.bulk = TRUE,
    group.by = "AST_status",
    ranges.title = region,
    ranges = li_enh,
    height.tracks = 20
  )
  
  annot = AnnotationPlot(object = g4,
                         region = region)
  ggsave(
    glue("{result_folder}{region}_cicero_covplots.pdf"),
    plot = p,
    width = 12,
    height = 8,
    device = "pdf"
  )
  ggsave(
    glue("{result_folder}{region}_annotplot.pdf"),
    plot = annot,
    width = 12,
    height = 8,
    device = "pdf"
  )
  return(print(p))
}

make_coverage_plot_all = function(region,
                                  distance = 5000,
                                  seurat_object = "../../results/Seurat/GFP_sorted_mousebrain/res0.8/outputs/Seurat_object.Rds",
                                  result_folder = "../../result_folder/") {
  
  g4 = readRDS(file = seurat_object)
  g4@meta.data = g4@meta.data %>%
    mutate(rownames_to_column(., var = "cell_id")) %>%
    mutate(AST_status = ifelse(cell_id %in% barcodes_scbridge, "AST", "non-AST"))
  
  pred_AST_links = ConnectionsToLinks(conns = AST, ccans = AST_ccan)
  Links(g4) = pred_AST_links
  AST_example = make_coverage_plot(region = region, seurat_object = g4)
  
  g4 = readRDS(file = seurat_object)
  g4@meta.data = g4@meta.data %>%
    mutate(rownames_to_column(., var = "cell_id")) %>%
    mutate(AST_status = ifelse(cell_id %in% barcodes_scbridge, "AST", "non-AST"))
  pred_nonAST_links = ConnectionsToLinks(conns = nonAST, ccans = nonAST_ccan)
  Links(g4) = pred_nonAST_links
  nonAST_example = make_coverage_plot(region = region, seurat_object = g4)
  
  plots = ggarrange(
    plotlist = list(nonAST_example, AST_example),
    ncol = 2,
    nrow = 1
  )
  
  ggsave(
    glue("{result_folder}{region}_cicero_arranged_covplots.pdf"),
    plot = plots,
    width = 12,
    height = 8,
    device = "pdf"
  )
  
}
