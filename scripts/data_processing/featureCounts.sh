# featureCount (FRIP)
PEAK=GSE173103_NAR_mESC_rep2_peaks.narrowPeak
BASE=GSE173103_NAR_mESC_rep2_peaks
BAM=../GSE173103_NAR_mESC_rep2.mm10.dedup.bam

# converting narrowPeak to SAF
awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print $4, $1, $2+1, $3, "."}' $PEAK > $BASE.saf

# FeatureCounts
featureCounts -p -a $BASE.saf -F SAF -o fCounts_$BASE.output.txt $BAM