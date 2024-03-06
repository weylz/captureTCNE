#!/bin/bash
# File    : call_peak.sh
# Time    : 2023/04/10 20:19:00
# Author  : Wenyong Zhu
# Version : 1.0.0
# Desc    : ...

[ -e GSE34448 ] || mkdir GSE34448
cd GSE34448
[ ! -e GSE34448_RAW.tar ] && echo 'ERROR! Please download GSE34448_RAW.tar from link (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE34448).'
tar xvf GSE34448_RAW.tar # download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE34448
mkdir cage_tags
mv GSM*_hg19_*.bigWig cage_tags/ && mv GSM*_hg19_*.bed.gz cage_tags/ && mv GSM*_hg19_*.bedRnaElements.gz cage_tags/
rm *.gtf.gz
# # A CAGE cluster is a region of overlapping tags with an assigned value that represents the expression level.
mkdir cage_peaks
num=0
total_num=`ls cage_tags/*.gz | wc -l`
for filename in `ls cage_tags/*.gz`
do
    fn=${filename#*/}
    fn=${fn%%.*}
    num=`expr $num + 1`
    echo [$num/$total_num] $fn
    less $filename | awk '{OFS=FS="\t";print $1,$2,$3,$1":"$2"-"$3,$7,$6}' > cage_peaks/$fn.bed
done
cd ..
cat ./GSE34448/cage_peaks/*bed | sortBed | mergeBed -i stdin -s -c 6 -o distinct | awk '{OFS=FS="\t";print $1,$2,$3,$1":"$2"-"$3,".",$4}' > ./peak_cage_GSE34448.bed

wget https://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/hg19.cage_peak_phase1and2combined_coord.bed.gz
cat ./Fantom5/hg19.cage_peak_phase1and2combined_coord.bed | awk '{OFS=FS="\t";print $1,$2,$3,$1":"$2"-"$3,".",$6}' > ./peak_cage_Fantom5.bed

cat peak_cage_GSE34448.bed peak_cage_Fantom5.bed | sortBed | mergeBed -i stdin -s -c 6 -o distinct | awk '{OFS=FS="\t";print $1,$2,$3,$1":"$2"-"$3,".",$4}' > peaks_cage_merged.bed
