suppressPackageStartupMessages({
  library("glue")
  library("ggplot2")
  library("data.table")
  library("tidyverse")
  library("data.table")
  library("ggalluvial")
  library("RColorBrewer")
})

# result folder
seurat_folder = "../../results/Seurat/"
scbr_folder = "../../results/scBridge/output/"

# Seurat integrations
unsorted_cl1 = readRDS("../../results/Seurat/unsorted_mousebrain/res0.1/outputs/Seurat_object_cl1.Rds")
unsorted_cl0 = readRDS("../../results/Seurat/unsorted_mousebrain/res0.1/outputs/Seurat_object_cl0.Rds")
sorted = readRDS("../../results/Seurat/GFP_sorted_mousebrain/res0.8/integration/outputs/G4_scRNA_integration.Rds")

unsorted_cl1_meta = unsorted_cl1@meta.data
unsorted_cl0_meta = unsorted_cl0@meta.data
sorted_meta = sorted@meta.data
uns_1_input = tibble("data" = "unsorted", seurat_cluster = paste0("unsorted ", unsorted_cl1_meta$seurat_clusters), label = "neuron", 
                     pred_score = unsorted_cl1_meta$predicted.cell_type.score) %>% 
  mutate(pred_score_cat = case_when(pred_score < 0.25 ~ "pred. score < 0.25", 
                                    0.25 < pred_score & pred_score < 0.50 ~  "0.25 < pred. score < 0.50",
                                    0.50 < pred_score & pred_score< 0.75 ~ "0.50 < pred. score < 0.75",
                                    0.75 < pred_score ~ "accurate"))
uns_0_input = tibble("data" = "unsorted", seurat_cluster = paste0("unsorted ", unsorted_cl0_meta$seurat_clusters), label = unsorted_cl0_meta$predicted.cell_type, 
                     pred_score = unsorted_cl0_meta$predicted.cell_type.score) %>% 
  mutate(pred_score_cat = case_when(pred_score < 0.25 ~ "pred. score < 0.25", 
                                    0.25 < pred_score & pred_score < 0.50 ~  "0.25 < pred. score < 0.50",
                                    0.50 < pred_score & pred_score< 0.75 ~ "0.50 < pred. score < 0.75",
                                    0.75 < pred_score ~ "accurate")) %>% 
  mutate(label = case_when(label == "Astrocytes" ~ "AST", 
                           label == "Oligodendrocytes" ~ "MOL", .default = label))
  

sort_input = tibble("data" = "GFP+", seurat_cluster = paste0("GFP+ ", sorted_meta$seurat_clusters), label = sorted_meta$pred_cell_type,
                    pred_score = sorted_meta$pred_max_score) %>% 
  mutate(label = case_when(label == "Astrocytes" ~ "AST", 
                           label == "Oligodendrocytes" ~ "MOL", .default = label)) %>% 
  mutate(pred_score_cat = case_when(pred_score < 0.25 ~ "pred. score < 0.25", 
                                    0.25 < pred_score & pred_score < 0.50 ~  "0.25 < pred. score < 0.50",
                                    0.50 < pred_score & pred_score< 0.75 ~ "0.50 < pred. score < 0.75",
                                    0.75 < pred_score ~ "accurate"))

inputs = rbind(sort_input, uns_1_input, uns_0_input)
inputs = inputs %>% group_by(data, seurat_cluster, label, pred_score_cat) %>% summarise(freq = n())

pred_score_cat_order = factor(inputs$pred_score_cat, levels = c("pred. score < 0.25", "0.25 < pred. score < 0.50", "0.50 < pred. score < 0.75", "accurate"))

# alluvial plot Seurat
ggplot(data = inputs,
       aes(axis1 = data, axis2 = seurat_cluster, axis3 = label,
           y = freq)) +
  scale_x_discrete(limits = c("sample", "Seurat cluster", "pred. label"), expand = c(.2, .05)) +
  xlab("") +
  geom_alluvium(aes(fill = pred_score_cat_order)) +
  scale_fill_viridis_d() +
  ylab("cell number") +
  labs(fill = "Seurat prediction") +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  )

ggsave(
  glue("{seurat_folder}Seurat_integrations-alluvial_plot.pdf"),
  width = 10,
  height = 10,
  dpi = 500,
)

# scBridge integrations
scbr_sorted_preds = fread("../../results/scBridge/output/Bartosovic_GFPsorted-scbridge_predictions.csv", header = TRUE) %>%
  mutate(scBridge_prediction = case_when(Prediction == "Astrocytes" ~ "AST", 
                           Prediction == "Oligodendrocytes" ~ "MOL", 
                           Prediction == "Novel (Most Unreliable)" ~ "unreliable", .default = Prediction)) %>% 
  dplyr::select(-Prediction)

scbr_sorted_meta = tibble(cell_id = rownames(sorted_meta), Seurat_cluster = sorted_meta$seurat_clusters, Seurat_prediction = sorted_meta$pred_cell_type) %>%
  mutate(Seurat_prediction = case_when(Seurat_prediction == "Astrocytes" ~ "AST", 
                           Seurat_prediction == "Oligodendrocytes" ~ "MOL", .default = Seurat_prediction))
  
scbr_sorted_meta = inner_join(scbr_sorted_preds, scbr_sorted_meta, by = c("V1" = "cell_id")) %>% 
  mutate(data = "GFP+") %>% 
  dplyr::select(data, Seurat_cluster, Seurat_prediction, scBridge_prediction)

inputs_sorted_scbr = scbr_sorted_meta %>% group_by(data, Seurat_cluster, Seurat_prediction, scBridge_prediction) %>% summarise(freq = n())

# alluvial plot Seurat vs. scBridge
ggplot(data = inputs_sorted_scbr,
       aes(
         axis1 = data,
         axis2 = Seurat_cluster,
         axis3 = Seurat_prediction,
         axis4 = scBridge_prediction,
         y = freq
       )) +
  scale_x_discrete(
    limits = c("sample", "Seurat cluster", "Seurat prediction", "scBridge prediction"),
    expand = c(.2, .05)
  ) +
  xlab("") +
  geom_alluvium(aes(fill = scBridge_prediction)) +
  scale_fill_manual(
    values = 
    c("AST" = "#fbb4ae",
    "COP-NFOL" = "#b3cde3",
    "OEC" = "#ccebc5",
    "MOL" = "#decbe4",
    "OPC" = "#fed9a6",
    "Pericytes" = "#ffffb3",
    "VEC" = "#e5d8bd",
    "VLMC" = "#fddaec",
    "unreliable" = "#d9d9d9"
  )) +
  ylab("cell number") +
  labs(fill = "Seurat prediction") +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text(size = 15, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  )

ggsave(
  glue("{scbr_folder}scBr_int_comp_w_Seurat-alluvial_plot.pdf"),
  width = 12,
  height = 10,
  dpi = 500,
)
