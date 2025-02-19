#!/bin/bash
#SBATCH -A naiss2025-22-67
#SBATCH -J sratools
#SBATCH -n 4
#SBATCH --mem=8GB
#SBATCH -t 24:00:00
#SBATCH -p shared

# work dir
cd /cfs/klemming/projects/snic/snic2020-6-3/SZABOLCS/scG4_study/data/bulk_CnT/G4_mESC

# MACS3 peak calling (python 3.9)
source /cfs/klemming/home/s/szhetey/.pyenv/versions/3.9.0/envs/MACS3_3.9/bin/activate

BAM=bowtie2_picard_output/G4_mESC_bulkCnT_R2.mm10.dedup.bam
BASE=G4_mESC_bulkCnT_R2

macs3 callpeak -t $BAM \
	--nomodel --shift -100 \
	--extsize 200 -q 0.05 -g 2652783500 \
	-n $BASE \
	--outdir macs3
	
