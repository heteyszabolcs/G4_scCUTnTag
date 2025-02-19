#!/bin/bash
#SBATCH -A naiss2024-22-108
#SBATCH -J pqs_deeptools
#SBATCH -n 4
#SBATCH --mem=8GB
#SBATCH -t 24:00:00
#SBATCH -p shared

cd /cfs/klemming/projects/snic/snic2020-6-3/SZABOLCS/scG4_study/git_repo/results/peak_overlaps

module load bioinfo-tools
module load ucsc-utilities/v421
module load bedtools
module load deepTools

# NAME="unsorted_common_peaks-cl1_cl0"

# awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$4}' ${NAME}.bed > ${NAME}.bedgraph
# sort -k1,1 -k2,2n ${NAME}.bedgraph > ${NAME}_sorted.bedgraph
# bedtools merge -i ${NAME}_sorted.bedgraph -c 1 -o count > counted
# awk '/\t1$/{print}' counted > filtered
# bedtools intersect -a ${NAME}_sorted.bedgraph -b filtered -wa > ${NAME}_sorted_filt.bedgraph
# bedGraphToBigWig ${NAME}_sorted_filt.bedgraph /cfs/klemming/projects/snic/snic2020-6-3/SZABOLCS/reference_data/mm10.chrom.size.tsv ${NAME}.bw

NAME="unsorted_PQS"

computeMatrix reference-point -o ${NAME}.mat.gz \
        -S unsorted_unique_peaks-cl0.bw unsorted_common_peaks-cl1_cl0.bw unsorted_unique_peaks-cl1.bw \
        -R unsorted_unique_peaks-cl0.bedgraph unsorted_cl1_cl0-intersection.bedgraph unsorted_unique_peaks-cl1.bedgraph \
        -b 3000 -a 3000 \
        --referencePoint center \
        --samplesLabel "PQS score" "PQS score" "PQS score"

gunzip ${NAME}.mat.gz
sed 's~nan~0~g' ${NAME}.mat > ${NAME}_.mat
gzip ${NAME}_.mat
mv ${NAME}_.mat.gz ${NAME}.mat.gz

plotHeatmap -m ${NAME}.mat.gz \
        -out ${NAME}.pdf \
        --refPointLabel "peak" \
        --heatmapHeight 14 \
        --yMin 0 0 0 \
        --yMax 40 40 40 \
        --zMin 0 0 0\
        -z "cluster 0" "both" "cluster 1" \
        --zMax 40 40 40 \
        --colorList "white, #756bb1" "white, #756bb1" "white, #756bb1" \
        --yAxisLabel "" \
        --xAxisLabel ""


rm *.bedgraph
rm counted
rm filtered
