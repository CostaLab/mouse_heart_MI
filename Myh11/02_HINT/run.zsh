#!/bin/sh
### Job name
#SBATCH -J Day10_2021-10-20_12:05:34.351535
#SBATCH -e logs/Day10_2021-10-20_12:05:34.351535.txt
#SBATCH -o logs/Day10_2021-10-20_12:05:34.351535.txt

### Time your job needs to execute, e. g. 15 min 30 sec
#SBATCH -t 100:00:00

### Memory your job needs per node, e. g. 1 GB
#SBATCH --mem=200G
#SBATCH --cpus-per-task=20

source ~/miniconda3/bin/activate
conda activate seurat3


aaaaa="Day10"
bbbbb="DS9"


cluster_file=./$aaaaa.txt
DS=$bbbbb

bam_file=/home/sz753404/data/project/Mouse_MI_ATAC/Myh11/data/${DS}/outs/possorted_bam.bam
output_location=../BAM

mkdir ../BAM

python split_bam_by_cluster.py $cluster_file $bam_file $output_location $aaaaa #forward
