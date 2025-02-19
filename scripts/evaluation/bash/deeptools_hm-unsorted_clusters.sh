#!/bin/bash
#SBATCH -A naiss2024-22-108
#SBATCH -J deeptools
#SBATCH -n 4
#SBATCH --mem=8GB
#SBATCH -t 24:00:00
#SBATCH -p shared

module load bioinfo-tools
module load deepTools

cd /cfs/klemming/projects/snic/snic2020-6-3/SZABOLCS/scG4_study/git_repo/results/peak_overlaps

computeMatrix reference-point -o unsorted_cl0_cl1.mat.gz \
        -S /cfs/klemming/projects/snic/snic2020-6-3/SZABOLCS/scG4_study/results/Seurat/final/unsorted_brain/res0.1/cluster_spec_bigwigs/0.bam_RPGC.bigwig \
		/cfs/klemming/projects/snic/snic2020-6-3/SZABOLCS/scG4_study/results/Seurat/final/unsorted_brain/res0.1/cluster_spec_bigwigs/1.bam_RPGC.bigwig \
		/cfs/klemming/projects/snic/snic2020-6-3/SZABOLCS/scG4_study/data/CellRanger/GFP_sorted/possorted_RPGC.bigwig \
        -R ../Seurat/unsorted_mousebrain/res0.1/cluster_spec_peaks/bedtools-unique_unsorted_cl0.bed \
		unsorted_cl1_cl0-intersection.bed \
		../Seurat/unsorted_mousebrain/res0.1/cluster_spec_peaks/bedtools-unique_unsorted_cl1.bed \
        -b 3000 -a 3000 \
        --referencePoint center \
        --samplesLabel "cluster 0" "cluster 1" "GFP+" \
		--skipZeros --missingDataAsZero
		
plotHeatmap -m unsorted_cl0_cl1.mat.gz \
        -out unsorted_cl0_cl1.pdf \
        --refPointLabel "peak" \
        --heatmapHeight 14 \
        --yMin 0 0 0 \
        --yMax 90 90 90 \
        --zMin 0 0 0 \
        -z "cluster 0" "both" "cluster 1" \
        --zMax 90 90 90 \
        --colorList "white, #31a354" "white, #31a354" "white, #31a354" \
        --yAxisLabel "" \
        --xAxisLabel ""