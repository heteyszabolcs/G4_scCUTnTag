#!/bin/bash
#SBATCH -A naiss2025-22-67
#SBATCH -J sratools
#SBATCH -n 4
#SBATCH --mem=8GB
#SBATCH -t 24:00:00
#SBATCH -p shared

# modules
module load bioinfo-tools
module load samtools
module load bowtie2
module load deepTools
module load picard

# go to dir of fastqs
cd /cfs/klemming/projects/snic/snic2020-6-3/SZABOLCS/scG4_study/data/GSE173103

# add bowtie2 index and fastq files
INDEX=/cfs/klemming/projects/snic/snic2020-6-3/SZABOLCS/reference_data/bowtie2_index/mm10/mm10
DATA_DIR=/cfs/klemming/projects/snic/snic2020-6-3/SZABOLCS/scG4_study/data/GSE173103/fastq
FASTQ1="GSE173103_NAR_mESC_rep1_R1.fastq.gz"
FASTQ2="GSE173103_NAR_mESC_rep1_R2.fastq.gz"

# workflow
# 1. alignment 2. bam converting 3. remove duplicates 4. bam coverage

base=$(basename $FASTQ1 _R1.fastq.gz)
echo $base
#echo $f
#bowtie2 -x $INDEX -p 2 --fast-local -1 $DATA_DIR/$FASTQ1 -2 $DATA_DIR/$FASTQ2 | samtools view -bS -F 4 -o $base.bam 
#echo "sorting $base.bam"
#samtools sort $base.bam -T $base -o $base.bam
echo "remove duplicates"
#samtools addreplacerg -r "@RG\tID:RG1\tSM:SampleName\tPL:Illumina\tLB:Library.fa" -o $base.picard.input.bam $base.bam
#java -jar $PICARD MarkDuplicates -REMOVE_DUPLICATES TRUE -I $base.picard.input.bam -O $base.mm10.dedup.bam -M $base.mm10.all.MarkDuplicates.txt
rm $base.picard.input.bam
#echo "indexing"
#samtools index $base.mm10.dedup.bam
echo "create RPGC bigwig"
#bamCoverage --bam $base.mm10.dedup.bam -o ${base}_RPGC.bigwig \
#   --binSize 10 \
#   --normalizeUsing RPGC \
#   --effectiveGenomeSize 2652783500

plotFingerprint -b bowtie2_picard_output/G4_mESC_bulkCnT_R1.mm10.dedup.bam \
	bowtie2_picard_output/G4_mESC_bulkCnT_R2.mm10.dedup.bam  \
	/cfs/klemming/projects/snic/snic2020-6-3/SZABOLCS/scG4_study/data/GSE173103/GSE173103_NAR_mESC_rep1.mm10.dedup.bam \
	/cfs/klemming/projects/snic/snic2020-6-3/SZABOLCS/scG4_study/data/GSE173103/GSE173103_NAR_mESC_rep2.mm10.dedup.bam \
	--labels "mESC - spin protocol (rep 1)" "mESC - spin protocol (rep 2)" "mESC - bead protocol, rep 1" "mESC - bead protocol, rep 2" \
	--minMappingQuality 20 --skipZeros \
	--region 19 --numberOfSamples 50000 \
	-T "Fingerprints of different bulk CUT&Tag protocols"  \
	--plotFile bulk_CnT-fingerprints.pdf 
