if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(
  "cicero",
  "monocle3",
  "Seurat",
  "Signac",
  "ggplot2",
  "data.table",
  "tidyverse",
  "glue",
  "GenomicRanges"
)

# if 'make_cicero_cds' function is not working:
# install monocle3 again with the line below BEFORE loading packages:
#devtools::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")
set.seed(42)

# cicero: https://www.bioconductor.org/packages/release/bioc/vignettes/cicero/inst/doc/website.html
# implementation: https://github.com/Castelo-Branco-lab/scCut-Tag_2020/blob/master/notebooks/H3K27ac/Cicero.Rmd

# result folder
result_folder = "../../results/astrocyte_exploration/"

# helper function
source("../../utils/cicero_helper.R")

# data
mm10.chromsize = fread("../../data/mm10.chrom.sizes.txt")
g4 = readRDS(file = "../../results/Seurat/GFP_sorted_mousebrain/res0.8/outputs/Seurat_object.Rds")
pred = fread("../../results/scBridge/output/scbridge_predictions.csv", header = TRUE)

# K27ac-ed mouse brain cCREs (Li et al.)
li_enh = readRDS("../../data/Li_et_al-mousebrain.union.cCRE_with_K27ac.Rds")
li_enh_bed = as_tibble(li_enh)
write_tsv(li_enh_bed, "../../data/Li_et_al-mousebrain_cCRE_with_K27ac.bed", col_names = FALSE)

# cicero workflow
input_cds = monocle3::new_cell_data_set(
  expression_data =  g4[['peaks']]@counts,
  cell_metadata = g4@meta.data,
  gene_metadata = data.frame(
    'gene_short_name' = rownames(g4[['peaks']]),
    row.names = rownames(g4[['peaks']])
  )
)
input_cds = detect_genes(input_cds)
input_cds = input_cds[Matrix::rowSums(exprs(input_cds)) != 0, ]

input_cds = detect_genes(input_cds)
input_cds = estimate_size_factors(input_cds)
input_cds = preprocess_cds(input_cds, method = "LSI")

input_cds = reduce_dimension(input_cds,
                             reduction_method = 'UMAP',
                             preprocess_method = "LSI")


umap_coords = reducedDims(input_cds)$UMAP
cicero_cds = make_cicero_cds(input_cds, reduced_coordinates = umap_coords)

# estimate the co-G4 sites in the genome in order to predict cis-regulatory interactions
conns = run_cicero(cicero_cds, mm10.chromsize)
save(conns, file = glue("{result_folder}cicero_GFPsorted.Rds"))

# finding cis-Coaccessibility networks (CCANS), coaccess_cutoff_override > 0.1 strict
CCAN_assigns = generate_ccans(conns, coaccess_cutoff_override = 0.1)
save(CCAN_assigns, file = glue("{result_folder}cicero_GFPsorted_coG4networks.Rds"))

# scBridge outputs
barcodes_scbridge = pred %>% filter(Prediction == "Astrocytes") %>% pull(V1)
other = setdiff(colnames(g4@assays$GA@counts), barcodes_scbridge)

g4@meta.data = g4@meta.data %>% 
  mutate(rownames_to_column(., var = "cell_id")) %>% 
  mutate(AST_status = ifelse(cell_id %in% barcodes_scbridge, "AST", "non-AST"))

g4_pred_ast = subset(g4, cells = barcodes_scbridge)
g4_pred_nonast = subset(g4, cells = other)

# cicero on AST predictions
pred_ast = cicero_function(seurat_object = g4_pred_ast, label = "predAST", 
                           output_path = result_folder)

# load(file = "../results/Seurat/callpeaks_GFPsorted/sorted_cicero-predAST.Rds")
# pred_ast_links = ConnectionsToLinks(conns = conns, ccans = pred_mol)
# Links(g4) = pred_mol_links

pred_nonast = cicero_function(seurat_object = g4_pred_nonast, label = "pred_nonAST",
                              output_path = result_folder)

##
load("../../results/astrocyte_exploration/cicero_GFPsorted-predAST.Rds")
load(
  "../../results/astrocyte_exploration/cicero_GFPsorted_coG4networks-predAST.Rds"
)
AST = conns
AST_ccan = CCAN_assigns
rm(conns)
load("../../results/astrocyte_exploration/cicero_GFPsorted-pred_nonAST.Rds")
load(
  "../../results/astrocyte_exploration/cicero_GFPsorted_coG4networks-pred_nonAST.Rds"
)
nonAST = conns
nonAST_ccan = CCAN_assigns
rm(conns)

# keep G4 peaks with coaccessibility above 0.2
coaccess_thr = 0.2

sign_coaccess_nonAST = nonAST %>%
  dplyr::filter(coaccess > coaccess_thr) %>%
  separate(Peak1,
           sep = "-",
           into = c("chr1", "start1", "end1")) %>%
  separate(Peak2,
           sep = "-",
           into = c("chr2", "start2", "end2")) %>%
  dplyr::select(chr1, start1, end2) %>%
  mutate(start1 = as.numeric(start1), end2 = as.numeric(end2)) %>%
  dplyr::filter(start1 < end2)

sign_coaccess_nonAST$type = "non-AST"
nonAST_peakset = GRanges(
  seqnames = sign_coaccess_nonAST$chr1,
  ranges = IRanges(
    start = sign_coaccess_nonAST$start1,
    end = sign_coaccess_nonAST$end2,
    names = sign_coaccess_nonAST$type,
  )
)

sign_coaccess_AST = AST %>%
  dplyr::filter(coaccess > coaccess_thr) %>%
  separate(Peak1,
           sep = "-",
           into = c("chr1", "start1", "end1")) %>%
  separate(Peak2,
           sep = "-",
           into = c("chr2", "start2", "end2")) %>%
  dplyr::select(chr1, start1, end2) %>%
  mutate(start1 = as.numeric(start1), end2 = as.numeric(end2)) %>%
  dplyr::filter(start1 < end2)

sign_coaccess_AST$type = "AST"
AST_peakset = GRanges(
  seqnames = sign_coaccess_AST$chr1,
  ranges = IRanges(
    start = sign_coaccess_AST$start1,
    end = sign_coaccess_AST$end2,
    names = sign_coaccess_AST$type,
  )
)

# keep AST specific regions
ol = findOverlaps(AST_peakset,
                  nonAST_peakset,
                  type = "any",
                  ignore.strand = FALSE)
AST_spec = as_tibble(AST_peakset[-queryHits(ol)])
AST_spec_regions = AST_spec %>% mutate(region = paste(seqnames, start, end, sep = "-")) %>%
  pull(region)

# find AST spec G4 regions over cCREs
ol_w_enh = findOverlaps(AST_peakset[-queryHits(ol)],
                        li_enh,
                        type = "any",
                        ignore.strand = FALSE)
AST_spec_enh_regions = as_tibble(AST_peakset[-queryHits(ol)][queryHits(ol_w_enh)])
AST_spec_enh_regions = distinct_all(AST_spec_enh_regions) %>%
  mutate(region = paste(seqnames, start, end, sep = "-")) %>%
  pull(region)

barcodes_scbridge = pred %>% filter(Prediction == "Astrocytes") %>% pull(V1)
other = setdiff(colnames(g4@assays$GA@counts), barcodes_scbridge)

g4@meta.data = g4@meta.data %>%
  mutate(rownames_to_column(., var = "cell_id")) %>%
  mutate(AST_status = ifelse(cell_id %in% barcodes_scbridge, "AST", "non-AST"))

# genome browser visualizations
selected_enhs = sample(AST_spec_enh_regions, 50)
lapply(selected_enhs, make_coverage_plot, result_folder = "../../results/cicero")

# example region
make_coverage_plot_all(region = "chr17-63492706-63500954", result_folder = "../../results/cicero/")
