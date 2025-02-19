suppressPackageStartupMessages({
  library("glue")
  library("tidyverse")
  library("data.table")
  library("wigglescout")
  library("ggrastr")
  library("pqsfinder")
  library("BSgenome")
  library("Biostrings")
  library("GenomicRanges")
  library("BSgenome.Mmusculus.UCSC.mm10")
})

# export folder
bed_folder = "../../results/peak_overlaps/"

# peak set
peaks = fread("../../results/peak_overlaps/unsorted_cl1_cl0-intersection.bed")
peaks = peaks %>% dplyr::filter(str_detect(V1, "GL|JH", negate = TRUE))

# reference genome (e.g. mm10)
mm10 = BSgenome.Mmusculus.UCSC.mm10
fasta = readDNAStringSet("C:/Szabolcs/Karolinska/Data/reference_data/genomes/mm10.fa")

# creating bed file containing PQS score (column 4th)
peaks$type = "peak"
gr = GRanges(
  seqnames = peaks$V1,
  ranges = IRanges(
    start = peaks$V2,
    end = peaks$V3,
    names = peaks$type,
  )
)

seq = Biostrings::getSeq(mm10, gr)

run_pqsf = function(x) {
  print(x)
  return(mean(score(pqsfinder(seq[x]$peak, min_score = 1))))
}
pqs_score = unlist(lapply(1:length(gr), run_pqsf))


bedg = peaks %>% dplyr::select(V1, V2, V3) %>% mutate(PQS_score = pqs_score) 
bedg[is.na(bedg$V4)]$V4 = 1
write_tsv(bedg, glue("{bed_folder}unsorted_common_peaks-cl1_cl0.bed"), col_names = FALSE)
