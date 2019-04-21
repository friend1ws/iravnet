#! /bin/bash
#! -S /bin/bash
#$ -cwd
#$ -e log/ -o log/

annot_utils boundary boundary.bed.gz --donor_size 3,6 --acceptor_size 6,1 --grc

annot_utils exon --grc exon.bed.gz

zcat db/TCGA_742_4.bed.gz | awk 'BEGIN{OFS="\t"} {print $1, $2-10, $2+10}' > control_4_broden.bed

bedtools subtract -a boundary.bed.gz -b control_4_broden.bed -A > boundary_proc.tmp.bed

bedtools subtract -a boundary_proc.tmp.bed -b exon.bed.gz -A -f 0.9 > boundary_proc.bed

rm -rf boundary.bed.gz*
rm -rf exon.bed.gz*
rm -rf control_4_broden.bed
rm -rf boundary_proc.tmp.bed

