if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(
  "tidyverse",
  "GenomicRanges",
  "org.Mm.eg.db",
  "org.Hs.eg.db",
  "TxDb.Mmusculus.UCSC.mm10.knownGene",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "ChIPseeker"
)

# hg38 human
hg38_annotation = function(regions,
                           seqname_col,
                           start_col,
                           end_col,
                           feature_1 = NULL,
                           feature_2 = NULL,
                           feature_3 = NULL) {
  
  seqnames = regions %>% pull(all_of(seqname_col))
  starts = regions %>% pull(all_of(start_col))
  ends = regions %>% pull(all_of(end_col))
  
  if(!is.null(feature_1)) {
    feature_1 = regions %>% pull(all_of(feature_1))
  }
  if(!is.null(feature_2)) {
    feature_2 = regions %>% pull(all_of(feature_1))
  }
  if(!is.null(feature_3)) {
    feature_3 = regions %>% pull(all_of(feature_1))
  }
  
  
  regions = GRanges(seqnames = seqnames,
                    ranges = IRanges(start = starts,
                                     end = ends,
                                     feature_1 = feature_1,
                                     feature_2 = feature_2))
  
  txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
  annot = annotatePeak(
    regions,
    tssRegion = c(-3000, 3000),
    TxDb = txdb,
    annoDb = "org.Hs.eg.db"
  )
  
  annot = as.data.frame(annot)
  return(annot)
  
}

# mm10 mouse
mm10_annotation = function(regions,
                           seqname_col,
                           start_col,
                           end_col,
                           feature_1 = NULL,
                           feature_2 = NULL,
                           feature_3 = NULL) {
  
  seqnames = regions %>% pull(all_of(seqname_col))
  starts = regions %>% pull(all_of(start_col))
  ends = regions %>% pull(all_of(end_col))
  
  if(!is.null(feature_1)) {
    feature_1 = regions %>% pull(all_of(feature_1))
  }
  if(!is.null(feature_2)) {
    feature_2 = regions %>% pull(all_of(feature_2))
  }
  if(!is.null(feature_3)) {
    feature_3 = regions %>% pull(all_of(feature_3))
  }
  
  
  regions = GRanges(seqnames = seqnames,
                    ranges = IRanges(start = starts,
                                     end = ends,
                                     feature_1 = feature_1,
                                     feature_2 = feature_2,
                                     feature_3 = feature_3))
  
  txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
  annot = annotatePeak(
    regions,
    tssRegion = c(-3000, 3000),
    TxDb = txdb,
    annoDb = "org.Mm.eg.db"
  )
  
  annot = as.data.frame(annot)
  return(annot)
  
}









