#! /bin/bash

CTYPE=$1

if [ ! -d output/${CTYPE} ]
then
    mkdir -p output/${CTYPE}
fi

ls /share/portal/TCGA/rna/${CTYPE}/hg19/alignment/star_2.5.2a/TCGA-*/*.bam > output/${CTYPE}/sample_list.txt

while read BAM
do
    BBAM=`basename ${BAM}`
    SAMPLE=${BBAM%.Aligned.sortedByCoord.out.bam}

    echo "qsub -l os6 irav_detect.tmp.sh ${BAM} output/${CTYPE}/${SAMPLE}.irav.result.txt"
    # qsub -l os6 irav_detect.tmp.sh ${BAM} output/${CTYPE}/${SAMPLE}.irav.result.txt
    bash irav_detect.tmp.sh ${BAM} output/${CTYPE}/${SAMPLE}.irav.result.txt

done < output/${CTYPE}/sample_list.txt

