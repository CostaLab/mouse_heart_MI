#!/bin/sh

#export RUBYLIB=$RUBYLIB:/home/rs619065/AMUSED:/home/rs619065/Ruby-DNA-Tools
#
ID=DS$1

case $1 in
  "1")
    sample="Healthy_Female"
    ;;

  "3")
    sample="Day3_Female"
    ;;

  "4")
    sample="Day3_Male"
    ;;

  "5")
    sample="Day10_Female"
    ;;

  "6")
    sample="Day10_Male"
    ;;


  "7")
    sample="Healthy_Male"
    ;;
esac



mkdir -p ../data/SplitBamFilesByClustering/BAM

cluster_file=../data/clusters/$sample.txt

bam_file=../../CountMatrix/$ID/outs/possorted_bam.bam
output_location=../data/SplitBamFilesByClustering/BAM
python split_bam_by_cluster.py $cluster_file $bam_file $output_location $sample

