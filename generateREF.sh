#!/bin/bash
# File    : generateREF.sh
# Time    : 2023/09/25 16:36:44
# Author  : Wenyong Zhu
# Version : 1.0.0
# Desc    : obtain full-length transcript and generate reference.


input_directory=$1
reference_directory=$2
suffix=$3

#Exit if file already exists
[ -e ../filtered_CNEs_BRCA.bed ] && echo "File filtered_CNEs_BRCA.bed already exists. Exit." && exit $?

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
unzip gencode.v19.annotation.gtf.gz
ref_gtf="$(pwd)/gencode.v19.annotation.gtf"

list=();for filename in $(ls $input_directory);do list+=($filename); done;num_sample=${#list[*]}

echo [ `date` ] "stringtie starting ..."
tmp_fifofile="/tmp/$$.fifo"; mkfifo $tmp_fifofile; exec 6<>$tmp_fifofile; rm $tmp_fifofile
thread_num=12; for ((i=0;i<${thread_num};i++));do echo; done >&6
num=0
for filename in ${list[*]}
do
	num=`expr $num + 1`
	tmp_list=`ls ./`
	if [[ "${tmp_list[@]}" =~ "$filename" ]]
	then
		echo [ `date` ] $num/$num_sample $filename "exists --> Skip."
	else
		read -u6
		{
			echo [ `date` ] $num/$num_sample $filename "stringtie - transcriptome assembly ..."
			stringtie -p 24 -G $ref_gtf -o $filename.gtf -l $filename $input_directory/$filename/$filename"_Aligned.sortedByCoord.out.sort.bam"
			echo [ `date` ] $num/$num_sample $filename "stringtie finished successfully."
			wait
			echo >&6
		} &
	fi
done
wait
exec 6>&-

echo [ `date` ] "stringtie merging ..."
for filename in $(ls $input_directory);do echo $filename.gtf >> list_gtf_assembly.txt; done
stringtie --merge -p 24 -G $ref_gtf -o stringtie_transcript_merged.gtf list_gtf_assembly.txt
cat stringtie_transcript_merged.gtf | awk '{FS=OFS="\t";if($3=="transcript") print $1,$4,$5,$1":"$4"-"$5,".",$7}' | grep "^chr" | sortBed > stringtie_transcript_merged.bed

echo [ `date` ] "generating full-length transcripts ..."
cat $ref_gtf | awk '{FS="\t|\"|\;";OFS="\t";if($16=="protein_coding"&&$3=="transcript") print $1,$4,$5,$1":"$4"-"$5,".",$7}' | sortBed > $reference_directory/gencode_protein_coding_transcript.bed
bedtools subtract -a stringtie_transcript_merged.bed -b $reference_directory/gencode_protein_coding_transcript.bed -f 1.0 -r | sortBed > stringtie_transcript_without_protein_coding.bed
echo [ `date` ] "Done."

echo [ `date` ] "generating reference ..."
cat $reference_directory/database/known_CNEs.bed | bedtools intersect -a stdin -b $reference_directory/cage/peaks_cage_merged.bed -u | bedtools intersect -a stdin -b stringtie_transcript_without_protein_coding.bed -f 1.0 -u | bedtools subtract -a stdin -b $reference_directory/blocklist/excluded_regions.bed -A | sortBed | mergeBed -i stdin | awk '{OFS="\t"}{if($3-$2>200) print $1,$2,$3,$1":"$2"-"$3}' > $reference_directory/filtered_CNEs_$suffix.bed
echo [ `date` ] "All Finished."
