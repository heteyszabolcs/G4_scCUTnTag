# create input
create_input = function(peaks) {
  annotations = mm10_annotation(regions = peaks, seqname_col = "V1", start_col = "V2", end_col = "V3", feature_1 = "V5", feature_2 = "V5", feature_3 = "V5")
  
  genes = annotations %>% dplyr::filter(abs(`distanceToTSS`) < 3000) %>%
    dplyr::filter(feature_1 > median(feature_1)) %>% 
    pull(SYMBOL) 
  
  input = rep(0, dim(annotations)[1])
  names(input) = annotations$SYMBOL
  input[names(input) %in% genes] = 1
  
  return(input)
}

# topGO analysis
create_go_matrix = function(genes, colname) {
  # find biological process ontology
  GOdata <- new(
    "topGOdata",
    ontology = "BP",
    # use biological process ontology
    allGenes = genes,
    geneSelectionFun = function(x)
      (x == 1),
    annot = annFUN.org,
    mapping = "org.Mm.eg.db",
    ID = "symbol"
  )
  
  # Fisher test
  resultFisher <-
    runTest(GOdata, algorithm = "elim", statistic = "fisher")
  out <-
    GenTable(GOdata,
             Fisher = resultFisher,
             topNodes = 10,
             numChar = 60)
  
  out = out %>% dplyr::select(Term, Fisher)
  colnames(out) = c("Term", colname)
  
  return(out)
  
}

# create_input function for FindMarkers outputs
create_input_fam = function(fc_table, background = peaks_1) {
  annotations = mm10_annotation(regions = background, seqname_col = "V1", start_col = "V2", end_col = "V3", feature_1 = "V5", feature_2 = "V5", feature_3 = "V5")

  background = annotations %>% dplyr::filter(abs(`distanceToTSS`) < 3000) %>%
    dplyr::filter(feature_1 > median(feature_1)) %>% 
    pull(SYMBOL) 
  
  genes = fc_table %>% pull(gene_symbol) 
  
  input = rep(0, length(background))
  names(input) = background
  input[names(input) %in% genes] = 1
  
  return(input)
}

# visualize enrichments by heatmap
generate_heatmap = function(matrix) {
  col_fun = colorRamp2(c(1e-10, 1e-12, 1e-16, 1e-20, 1e-30), c("grey", "#ece7f2", "#e5f5f9", "#99d8c9", "#2ca25f"))
  hm = Heatmap(
    matrix,
    column_title = "",
    row_title = "",
    name = "Fisher test",
    col = col_fun,
    heatmap_legend_param = list(
      title = "Fisher test",
      at = c(1, -0.05, -0.001),
      labels = c("not enriched", "p < 1e-20", "p < 1e-10"),
      legend_height = unit(2, "cm")
    ),
    rect_gp = gpar(col = "black", lwd = 0.1),
    show_column_dend = FALSE,
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    show_row_dend = FALSE,
    heatmap_width = unit(7, "cm"),
    heatmap_height = unit(10, "cm"),
    row_names_gp = gpar(fontsize = 5),
    column_names_gp = gpar(fontsize = 6),
    column_names_rot = 90
  )
  hm
}

generate_heatmap_fm = function(matrix) {
  col_fun = colorRamp2(c(1, 0.1, 0.05, 1e-2, 1e-3), c("grey", "#ece7f2", "#e5f5f9", "#99d8c9", "#2ca25f"))
  hm = Heatmap(
    matrix,
    column_title = "",
    row_title = "",
    name = "Fisher test",
    col = col_fun,
    heatmap_legend_param = list(
      title = "Fisher test",
      at = c(1, -0.05, -0.001),
      labels = c("not enriched", "p < 1e-3", "p < 1e-2"),
      legend_height = unit(2, "cm")
    ),
    rect_gp = gpar(col = "black", lwd = 0.1),
    show_column_dend = FALSE,
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    show_row_dend = FALSE,
    heatmap_width = unit(7, "cm"),
    heatmap_height = unit(10, "cm"),
    row_names_gp = gpar(fontsize = 5),
    column_names_gp = gpar(fontsize = 6),
    column_names_rot = 90
  )
  hm
}