# packages
print("Load R packages")
if (!require("pacman"))
  install.packages("pacman")
pacman::p_load("tidyverse", "data.table", "ggplot2", "glue", "wigglescout")

# featureCounts (FRIP)
samples = c("mESC bead CnT, rep1", "mESC bead CnT, rep2", "mESC spin CnT, rep1", "mESC spin CnT, rep2")
frips = c(6.4, 6.4, 22.6, 23.4) # FRIP scores coming from the featureCounts runs
input = tibble(sample = samples, frip = frips)

frip_bar = ggplot(data = input, aes(x = sample, y = frip, fill = sample)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c("#fcae91", "#fcae91", "#de2d26", "#de2d26")) +
  ylim(0, 50) +
  labs(
    title = "FRIP scores of bulk CUT&Tag protocols",
    x = "",
    y = "Fraction of reads in peaks"
  ) +
  guides(fill="none") +
  theme_classic() +
  theme(
    text = element_text(size = 9),
    plot.title = element_text(size = 13),
    axis.title.y = element_text(size = 13, color = "black"),
    axis.text.x = element_text(size = 10, color = "black", angle=90, vjust=0.5),
    axis.text.y = element_text(size = 10, color = "black")
  ) 
frip_bar

ggsave(
  glue("../../results/bulk_CnT_qc/mESC_bulk_CnT_FRIPs.pdf"),
  plot = frip_bar,
  width = 5,
  height = 5
)

# density scatter (wigglescout)
bws = c(
  "../../data/bulk_CUTnTag/bulk_G4_CnT_mESC_rep1_RPGC.bigwig",
  "../../data/bulk_CUTnTag/bulk_G4_CnT_mESC_rep2_RPGC.bigwig",
  "../../data/bulk_CUTnTag/bulk_G4_CnT_GSE173103_NAR_mESC_rep1_RPGC.bigwig",
  "../../data/bulk_CUTnTag/bulk_G4_CnT_GSE173103_NAR_mESC_rep2_RPGC.bigwig"
)

dens_sc = function(bigwig1, bigwig2) {
  bins_10kb = bw_bins(
    c(bigwig1, bigwig2),
    bin_size = 10000,
    genome = "mm10"
  )
  bins_10kb = as.data.frame(bins_10kb)
  bins_10kb = na.omit(bins_10kb)
  
  rep1 = colnames(bins_10kb)[6]
  rep2 = colnames(bins_10kb)[7]
  label1 = strsplit(rep1, "_RPGC")[[1]][1]
  label2 = strsplit(rep2, "_RPGC")[[1]][1]
  
  bins_10kb = bins_10kb %>% mutate(x = log2(.[[rep1]]), y = log2(.[[rep2]]))
  
  p = ggplot(bins_10kb, aes(x = x, y = y)) +
    geom_hex(bins = 1000) +
    scale_fill_continuous(type = "viridis") +
    xlab(label1) +
    ylab(label2) +
    xlim(-10, 10) +
    ylim(-10, 10) +
    theme_classic() +
    theme(
      text = element_text(size = 20),
      plot.title = element_text(size = 15),
      axis.text.x = element_text(size = 13, color = "black")
    )
  print(p)
  
  ggsave(
    glue("../../results/bulk_CnT_qc/{label1}-{label2}.pdf"),
    plot = last_plot(),
    width = 10,
    height = 10
  )
  
  return(p)
}

# run through
l = list()
for (i in seq(1, length(bws))) {
  for (j in seq(1, length(bws))) {
    p = dens_sc(bigwig1 = bws[i], bigwig2 = bws[j])
    l[[paste0(as.character(i), "_", as.character(j))]] = p
  }
}
