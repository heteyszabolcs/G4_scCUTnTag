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

# source annotation script for mm10
source("../../utils/annotation.R")

# export folder
result_folder = "../../results/Seurat/unsorted_mousebrain/"

## unsorted G4 scCUT&Tag
# unique peaks
unique_cl0 = fread(
  "../../results/Seurat/unsorted_mousebrain/res0.1/cluster_spec_peaks/bedtools-unique_unsorted_cl0.bed"
)
unique_cl1 = fread(
  "../../results/Seurat/unsorted_mousebrain/res0.1/cluster_spec_peaks/bedtools-unique_unsorted_cl1.bed"
)

# enrichR analysis on the unique unsorted cluster 1 peaks
unique_cl1 = unique_cl1 %>% dplyr::filter(V7 > quantile(V7, .90))
unique_cl1 = mm10_annotation(
  regions = unique_cl1,
  seqname_col = "V1",
  start_col = "V2",
  end_col = "V3"
)
unique_cl1_input = unique_cl1 %>%
  dplyr::filter(abs(distanceToTSS) < 3000)

dbs = c(
  "GO_Molecular_Function_2023",
  "GO_Cellular_Component_2023",
  "GO_Biological_Process_2023"
)
unique_cl1 = enrichr(unique_cl1_input$SYMBOL, dbs)
unique_cl1_biol = unique_cl1[["GO_Biological_Process_2023"]]
unique_cl1_biol = unique_cl1_biol %>% dplyr::filter(0.05 > P.value)

# horizontal bar for biol. processes
unique_cl1_bars = ggplot(unique_cl1_biol[1:30, ], aes(x = reorder(Term, -log10(P.value)), y = -log10(P.value))) +
  geom_bar(stat = 'identity',
           fill = "#9ecae1",
           color = "black") +
  coord_flip() +
  labs(title = "unique peaks of unsorted cluster 1", y = "-log10(p-value)", x = "GO Biologocial process (2023)") +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 14, color = "black")
  )
unique_cl1_bars
unique_cl1_terms = tibble("Term" = unique_cl1_biol$Term,
                          "p-value" = unique_cl1_biol$P.value,
                          type = "unique_cl1_terms")

# enrichR analysis on the unique unsorted cluster 0 peaks

unique_cl0 = unique_cl0 %>% dplyr::filter(V7 > quantile(V7, .90))
unique_cl0 = mm10_annotation(
  regions = unique_cl0,
  seqname_col = "V1",
  start_col = "V2",
  end_col = "V3"
)
unique_cl0_input = unique_cl0 %>%
  dplyr::filter(abs(distanceToTSS) < 3000)

dbs = c(
  "GO_Molecular_Function_2023",
  "GO_Cellular_Component_2023",
  "GO_Biological_Process_2023"
)
unique_cl0 = enrichr(unique_cl0_input$SYMBOL, dbs)
unique_cl0_biol = unique_cl0[["GO_Biological_Process_2023"]]
unique_cl0_biol = unique_cl0_biol %>% dplyr::filter(0.05 > P.value)

# horizontal bar for biol. processes
unique_cl0_bars = ggplot(unique_cl0_biol[1:30, ], aes(x = reorder(Term, -log10(P.value)), y = -log10(P.value))) +
  geom_bar(stat = 'identity',
           fill = "#9ecae1",
           color = "black") +
  coord_flip() +
  labs(title = "unique peaks of unsorted cluster 1", y = "-log10(p-value)", x = "GO Biologocial process (2023)") +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 14, color = "black")
  )
unique_cl0_bars
unique_cl0_terms = tibble("Term" = unique_cl0_biol$Term,
                          "p-value" = unique_cl0_biol$P.value,
                          type = "unique_cl0_terms")

# keep CNS related GO terms
terms = rbind(unique_cl0_terms, unique_cl1_terms) %>%
  mutate(Term = str_to_lower(Term)) %>%
  dplyr::filter(
    str_detect(
      Term,
      "nervous|neural|oligodendrocyte|microglia|neuron|axon|synaptic|synapsis|myelination|glial|glia"
    )
  )

terms = terms %>% pivot_wider(., names_from = "type", values_from = "p-value")
terms[is.na(terms)] = 1
terms = as.data.frame(terms)
rownames(terms) = terms$Term
terms = terms %>% dplyr::select(-Term)

# visualize p-values
# heatmap of p-values
col_fun = colorRamp2(c(1, 0.05, 0.025, 0.001, 0.002),
                     c("grey", "#fee0d2", "#fc9272", "#de2d26", "#99000d"))
hm = Heatmap(
  terms,
  column_title = "enrichR",
  row_title = "",
  name = "p-value",
  col = col_fun,
  heatmap_legend_param = list(
    title = "p-value",
    at = c(0.10, 0.05, 0.001),
    labels = c("not enriched", "p < 0.05", "p < 0.001"),
    legend_height = unit(2, "cm")
  ),
  rect_gp = gpar(col = "black", lwd = 0.1),
  show_column_dend = FALSE,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  heatmap_width = unit(7, "cm"),
  heatmap_height = unit(10, "cm"),
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 6),
  column_names_rot = 90
)
hm

# export heatmap
pdf(
  file = glue("{result_folder}enrichR-unsorted_cl_enrichr.pdf"),
  width = 5,
  height = 5
)
print(hm)
dev.off()
