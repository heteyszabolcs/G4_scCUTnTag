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
  
  # finding cis-Coaccessibility networks (CCANS), coaccess_cutoff_override > 0.5 strict
  CCAN_assigns = generate_ccans(conns, coaccess_cutoff_override = 0.1)
  save(CCAN_assigns, file = glue("{output_path}cicero_GFPsorted_coG4networks-{label}.Rds"))
  
  return(CCAN_assigns)
}