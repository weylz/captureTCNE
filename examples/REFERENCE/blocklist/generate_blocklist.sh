#!/bin/bash
# File    : generate_blocklist.sh
# Time    : 2023/09/25 10:07:15
# Author  : Wenyong Zhu
# Version : 1.0.0
# Desc    : ...

cat gencode_v19/gencode.v19.*.bed RNAcentral/RNAcentral_*.bed ucsc/ucsc_*.bed | cut -f1-6 | sortBed | mergeBed -i stdin -s -c 6 -o distinct | awk '{OFS=FS="\t";print $1,$2,$3,$1":"$2"-"$3,".",$4}' > excluded_regions.bed
