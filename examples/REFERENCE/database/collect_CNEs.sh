#!/bin/bash
# File    : collect_CNEs.sh
# Time    : 2023/09/25 10:01:43
# Author  : Wenyong Zhu
# Version : 1.0.0
# Desc    : ...

[ -e Ancora ] || mkdir Ancora
cd Ancora
wget http://ancora.genereg.net/downloads/hg38/vs_mouse/HCNE_hg38_mm10_80pc_50col.bed.gz && gunzip HCNE_hg38_mm10_80pc_50col.bed.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
sed '1d' HCNE_hg38_mm10_80pc_50col.bed | liftOver stdin hg38ToHg19.over.chain.gz HCNE_hg19_mm10_80pc_50col.bed HCNE_hg19_mm10_80pc_50col.unMapped
rm HCNE_hg19_mm10_80pc_50col.unMapped
wget http://ancora.genereg.net/downloads/hg19/vs_rat/HCNE_hg19_rn4_80pc_50col.bed.gz
wget http://ancora.genereg.net/downloads/hg19/vs_rat/HCNE_hg19_rn5_80pc_50col.bed.gz
gunzip HCNE_hg19_rn4_80pc_50col.bed.gz HCNE_hg19_rn5_80pc_50col.bed.gz
sed '1d' HCNE_hg19_rn4_80pc_50col.bed > HCNE_hg19_rnx_80pc_50col.bed && sed '1d' HCNE_hg19_rn5_80pc_50col.bed >> HCNE_hg19_rnx_80pc_50col.bed
rm HCNE_hg19_rn4_80pc_50col.bed HCNE_hg19_rn5_80pc_50col.bed
intersectBed -a HCNE_hg19_mm10_80pc_50col.bed -b HCNE_hg19_rnx_80pc_50col.bed | cut -f1-3 | sortBed | mergeBed | awk '{OFS="\t"}{if($3-$2>=200) print $1,$2,$3,$1"_"$2"_"$3}' | sortBed > HCNE.Ancora.bed
cd ..

[ -e UCNEbase ] || mkdir UCNEbase
cd UCNEbase
wget https://epd.expasy.org/ucnebase/data/download/ucnes/hg19_UCNE_coord.bed
cat hg19_UCNE_coord.bed | cut -f1-3 | sortBed | mergeBed -i stdin > UCNE.UCNEbase.bed
cd ..

cat Ancora/HCNE.Ancora.bed UCNEbase/UCNE.UCNEbase.bed | sortBed | mergeBed -i stdin | awk '{FS=OFS="\t";print $1,$2,$3,$1":"$2"-"$3}' > known_CNEs.bed
