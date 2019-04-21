#! /bin/bash
#! -S /bin/bash
#$ -cwd
#$ -e log/ -o log/

annot_utils boundary boundary.bed.gz --donor_size 3,6 --acceptor_size 6,1 --grc

annot_utils exon --grc exon.bed.gz

zcat control_4.bed.gz | awk 'BEGIN{OFS="\t"} {print $1, $2-10, $2+10}' > control_4_broden.bed

bedtools subtract -a boundary.bed.gz -b contron_4_broden.bed > boundary_proc.tmp.bed

bedtools subtract -a boundary_proc.tmp.bed -b exon.bed.gz > boundary_proc.bed

rm -rf boundary_proc.tmp.bed

