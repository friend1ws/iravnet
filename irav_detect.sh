#! /bin/bash
#! -S /bin/bash
#$ -cwd
#$ -e log/ -o log/

export PYTHONHOME=/usr/local/package/python/current2.7
export PYTHONPATH=~/.local/lib/python2.7/site-packages
export PATH=${PYTHONHOME}/bin:~/.local/bin:$PATH
export PATH=${PYTHONHOME}/bin:$PATH
export LD_LIBRARY_PATH=${PYTHONHOME}/lib:${LD_LIBRARY_PATH}

export LANG=en_US.UTF-8

export PATH=/home/yshira/bin/htslib-1.7:${PATH}
export PATH=/home/yshira/bin/samtools-1.9:${PATH}
export PATH=/home/yshira/bin/bedtools2/bin:$PATH
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/yshira/bin/Complete-Striped-Smith-Waterman-Library/src

INPUT=$1
OUTPUT=$2

REF=/home/w3varann/database/GRCh37/GRCh37.fa

samtools mpileup ${INPUT} -l boundary_proc.bed -f ${REF} -O > ${OUTPUT}.tmp1.txt

python proc_variant.py ${OUTPUT}.tmp1.txt > ${OUTPUT}.tmp2.txt

intron_retention_utils allele_count ${INPUT} ${OUTPUT}.tmp2.txt ${OUTPUT}.tmp3.txt --reference ${REF}

python filter_irav.py ${OUTPUT}.tmp3.txt > ${OUTPUT}

# rm -rf ${OUTPUT}.tmp1.txt
# rm -rf ${OUTPUT}.tmp2.txt
# rm -rf ${OUTPUT}.tmp3.txt

