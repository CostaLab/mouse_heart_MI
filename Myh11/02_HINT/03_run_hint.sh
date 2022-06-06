#!/bin/bash
#
#### Job name
#SBATCH -J myo_2021-10-20_13:43:53.434189
#SBATCH -e logs/myo_2021-10-20_13:43:53.434189.txt
#SBATCH -o logs/myo_2021-10-20_13:43:53.434189.txt
#SBATCH -t 120:00:00
#SBATCH --mem=100G --cpus-per-task=40

source ~/miniconda3/bin/activate
conda activate seurat3


PDIR=`pwd`

bam_loc=$PDIR/DiffFootprinting/BAM
peaks_loc=$PDIR/DiffFootprinting/Peaks
footprint_loc=$PDIR/DiffFootprinting/Footprints
motifmatching_loc=$PDIR/DiffFootprinting/MotifMatching
diff_loc=$PDIR/DiffFootprinting/Diff
bigwig_loc=$PDIR/DiffFootprinting/Diff

mkdir -p ${peaks_loc}
mkdir -p ${footprint_loc}
mkdir -p ${bigwig_loc}
mkdir -p ${motifmatching_loc}
mkdir -p ${bam_loc}
mkdir -p ${diff_loc}

echo ${footprint_loc}

cp  ../BAM/*Myofibroblasts.bam $PDIR/DiffFootprinting/BAM



footprinting(){
    	filename=$1
    	samtools index $1
    	macs2 callpeak -t $filename -n ${filename%.bam} --outdir $2 -g mm --nomodel -f BAMPE -q 0.01 --keep-dup all
    	sed -i '/random\|chrUn\|chrM/d' $2/${filename%.bam}_peaks.narrowPeak
    	sed -i '/random\|chrUn\|chrM/d' $2/${filename%.bam}_summits.bed
	rgt-hint footprinting --organism=mm10 --atac-seq --paired-end --output-location=$4 --output-prefix=${filename%.bam} $1 $2/${filename%.bam}_peaks.narrowPeak
	rgt-motifanalysis matching --organism=mm10 --output-location=$5 --input-files $4/${filename%.bam}.bed
}

cd ${bam_loc}
for filename in *.bam;
do
	footprinting $filename ${peaks_loc} ${bigwig_loc} ${footprint_loc} ${motifmatching_loc} &
done

wait

cd ${diff_loc}
rgt-hint differential --organism=mm10 --bc --nc 40 \
--mpbs-files=${motifmatching_loc}/Sham_Myofibroblasts_mpbs.bed,\
${motifmatching_loc}/Day3_Myofibroblasts_mpbs.bed,\
${motifmatching_loc}/Day10_Myofibroblasts_mpbs.bed \
--reads-files=${bam_loc}/Sham_Myofibroblasts.bam,\
${bam_loc}/Day3_Myofibroblasts.bam,\
${bam_loc}/Day10_Myofibroblasts.bam \
--conditions=Sham,Day3,Day10 \
--output-location=${diff_loc} \
--output-prefix=All

