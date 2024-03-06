#!/bin/bash
# File    : captureTCNE.sh
# Time    : 2023/01/01 00:00:00
# Author  : Wenyong Zhu
# Version : 1.0.0
# Desc    : To capture the transcribed conserved noncoding elements in sample group according to the defined criteria.


# Requirement
#     - Python 3.6 or later;
#     - R 4.0 or later (library fitdistrplus 1.1-8);
#     - BEDTools suite v2.30.0;
#     - Jim Kent's executable programms: http://hgdownload.cse.ucsc.edu/admin/exe/.

# $@-List of all parameters
Options=$@
# $#-The number of arguments added to the Shell
Optnum=$#

# $0-The file name of the Shell itself ($0 $Options)
function usage(){
cat << EOF

Usage:    ${0##*\/} -g <...> -r <...> [-m <...>] [-c <...>] [-v] [-h]

    Note: -g and -r are required.

Requirement:
    - Python 3.6 or later;
    - R 4.0 or later (R library - fitdistrplus);
    - BEDTools suite;
    - stringTie;
    - Jim Kent's executable programms: http://hgdownload.cse.ucsc.edu/admin/exe/.

Options:
    -h    Print this help menu.
    -g    Label for sample group. The name of the sample group label (the name of the folder where the bigwig files are placed).
    -r    Absolute Path of reference element file. The bed file without a header requires four columns for chrom, chromStart, chromEnd, and name, respectively.
    -m    Minimum number (included) of transcriptional signals in the sample group, default: 1.
    -c    Excluding samples of sample group with TPM values in the lower ... of the normal distribution. Range: [0,1), default: 0.05.
    -v    The version of this tool.

EOF
}

function scriptInfo(){
    if [[ $1 == "errorGL" ]]; then
        echo "--------------------------------------------------"
        echo ">>>>>>>>>>>>> Invalid Syntax [-g ..] <<<<<<<<<<<<<"
        echo "--------------------------------------------------"
    elif [[ $1 == "errorRD" ]]; then
        echo "--------------------------------------------------"
        echo ">>>>>>>>>>>>> Invalid Syntax [-r ..] <<<<<<<<<<<<<"
        echo "--------------------------------------------------"
    elif [[ $1 == "errorNM" ]]; then
        echo "--------------------------------------------------"
        echo ">>>>>>>>>>>>> Invalid Syntax [-m ..] <<<<<<<<<<<<<"
        echo "--------------------------------------------------"
    elif [[ $1 == "errorCO" ]]; then
        echo "--------------------------------------------------"
        echo ">>>>>>>>>>>>> Invalid Syntax [-c ..] <<<<<<<<<<<<<"
        echo "--------------------------------------------------"
    elif [[ $1 == "errorExcept" ]]; then
        echo "--------------------------------------------------"
        echo ">>>>>>>>>>>>> Invalid Syntax [-* ..] <<<<<<<<<<<<<"
        echo "--------------------------------------------------"
    elif [[ $1 == "start" ]]; then
        echo "--------------------------------------------------"
        echo ">>>>>>>>>>>>>>>>>> Starting ... <<<<<<<<<<<<<<<<<<"
        echo "--------------------------------------------------"
    elif [[ $1 == "finish" ]]; then
        echo "--------------------------------------------------"
        echo ">>>>>>>>>>>>>>>>>> ALL FINISHED <<<<<<<<<<<<<<<<<<"
        echo "--------------------------------------------------"
    elif [[ $1 == "version" ]]; then
        echo -e "captureTCNE v1.0.0\n"
    fi
}

# $?-Closing code of the last command to run (return value)
if [ $# = 0 ]; then scriptInfo errorExcept && usage && exit $?; fi

# ;;-Use the options in case and play the role of Terminator
while getopts ":hg:r:m:c:v" OPTION
do
    case "$OPTION" in
        h ) usage && exit $?                   ;;
        g ) group_label="$OPTARG"              ;;
        r ) ref_data="$OPTARG"                 ;;
        m ) num_min="$OPTARG"                  ;;
        c ) cut_off="$OPTARG"                  ;;
        v ) scriptInfo version && exit $?      ;;
        * ) scriptInfo errorExcept && exit $?  ;;
    esac
done

# The name of the sample group label (the name of the folder where the bigwig files are placed).
[ -z "$group_label" ] && [ ! -d "$group_label"  ] && scriptInfo errorGL && exit $?
# The bed file without a header requires four columns for chrom, chromStart, chromEnd, and name, respectively.
[ -z "$ref_data" ] && [ ! -e "$ref_data" ] && scriptInfo errorRE && exit $?
# The minimum number (included) of transcriptional signals in the sample group (default: 1).
[ -z "$num_min" ] && num_min=1
# Excluding samples of sample group with TPM values in the lower ... of the normal distribution. Range: (0,1), default: 0.05.
[ -z "$cut_off" ] && cut_off=0.05
if [[ `echo "$cut_off>=1" | bc` -eq 1 ]] || [[ `echo "$cut_off<0" | bc` -eq 1 ]]; then scriptInfo errorCO && exit $?; fi

#Exit if file already exists
[ -e "TCNE."$group_label".bed" ] && echo "File TCNE."$group_label".bed already exists. Exit." && exit $?

echo -e "\nScript: ${0##*\/} $@ \n" && scriptInfo start

cd $group_label

echo "["`date`"] scanning files in sample group ..."
list_sample=()
for filename in $(ls ./)
do
    if [[ $filename == *.bw ]]
    then
        list_sample+=(${filename%.*})
    fi
done
num_sample=${#list_sample[*]}

[ $num_min -gt $num_sample ] && scriptInfo errorNM && exit $?

[ $num_sample == 0 ] && echo -e "No files in sample group.\nInterrupt!" && scriptInfo finish && exit $?

[ -e transcriptional.signal.cutoff.value.txt ] && : > transcriptional.signal.cutoff.value.txt
echo ["`date`"] "merging transcriptional signal of sample group ..."
[ -e bwfiles ] || mkdir bwfiles && mv *.bw bwfiles
bigWigMerge bwfiles/*.bw $group_label.bedGraph

echo ["`date`"] "measuring transcriptional signal distribution for sample group ..."
fitTSD.R $group_label.bedGraph $cut_off

echo ["`date`"] "converting file format in sample group ..."
[ -e bgfiles ] || mkdir bgfiles
rm -rf temp_files && mkdir temp_files
count=0
for filename in ${list_sample[*]}
do
    count=`expr $count + 1`
    echo -e "\t"$count/$num_sample "-->" $filename "..."
    bigWigToBedGraph bwfiles/$filename.bw temp_files/$filename.bedGraph
    sed -i '/^GL\|^JH\|^KB\|^KE/d' temp_files/$filename.bedGraph
    if [ `grep -c "chr" temp_files/$filename.bedGraph` -eq 0 ]
    then
        cat temp_files/$filename.bedGraph | awk '{print "chr"$0}' | sortBed -i > bgfiles/$filename.bedGraph
    else
        cat temp_files/$filename.bedGraph | sortBed -i > bgfiles/$filename.bedGraph
    fi
done
rm -rf temp_files

echo ["`date`"] "calculating transcriptional signal of elements in sample group ..."
[ -e bedfiles ] || mkdir bedfiles
count=0
for filename in ${list_sample[*]}
do
    count=`expr $count + 1`
    echo -e "\t"$count/$num_sample "-->" $filename "..."
    mapBed -a $ref_data -b bgfiles/$filename.bedGraph -c 4 -o max > bedfiles/$filename.re.bed
done

echo ["`date`"] "generating transcriptional signal matrix of elements in each sample ..."
[ -e tsfiles ] || mkdir tsfiles
count=0
for filename in ${list_sample[*]}
do
    count=`expr $count + 1`
    echo -e "\t"$count/$num_sample "-->" $filename "..."
    cut -f5 bedfiles/$filename.re.bed | sed "1i\\$filename" > tsfiles/$filename.tmp
done

echo ["`date`"] "integrating transcriptional signal of sample group ..."
cut -f1-4 bedfiles/${list_sample[0]}.re.bed | sed "1i\chrom\\tchromStart\\tchromEnd\\tname" > tmp_header.tmp
paste tmp_header.tmp tsfiles/*.tmp > transcriptional_signal.matrix
mv tmp_header.tmp tsfiles

echo ["`date`"] "filtering eligible elements from transcriptional signal matrix of sample group ..."
filterMatrix.py $num_sample $num_min
cat TCNE.over"$num_min".matrix | awk '{OFS="\t"} BEGIN{sum = 0} {for(i = 5; i <= NF; i++) {sum += $i} {print $1,$2,$3,$1"_"$2"_"$3, sum/(NF-4), "."; sum=0}}' | sed '1d' > TCNE."$group_label".bed
echo ["`date`"] "Please check the file: TCNE."$group_label".bed."

scriptInfo finish
