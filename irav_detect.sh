#! /bin/bash
#! -S /bin/bash
#$ -cwd
#$ -e log/ -o log/

INPUT=$1
OUTPUT=$2

REF=/home/w3varann/database/GRCh37/GRCh37.fa

samtools mpileup ${INPUT} -l boundary_proc.bed -f ${REF} -O > ${OUTPUT}.tmp1.txt

python proc_variant.py ${OUTPUT}.tmp1.txt > ${OUTPUT}.tmp2.txt

intron_retention_utils allele_count ${INPUT} ${OUTPUT}.tmp2.txt ${OUTPUT}.tmp3.txt --reference ${REF}

python filter_irav.py ${OUTPUT}.tmp3.txt > ${OUTPUT}

rm -rf ${OUTPUT}.tmp1.txt
rm -rf ${OUTPUT}.tmp2.txt
rm -rf ${OUTPUT}.tmp3.txt

