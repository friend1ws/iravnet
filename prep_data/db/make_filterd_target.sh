#! /usr/bin/env bash

# hg19 grc
annot_utils boundary boundary.bed.gz --donor_size 3,6 --acceptor_size 6,1 --grc

# annot_utils exon exon.bed.gz --grc

zcat TCGA_742_8.bed.gz | awk '$1 ~ /^[0-9XY]+$/ {OFS="\t"; print $1, $2-10, $2+10}' > control_8_broden.bed

# bedtools subtract -a boundary.bed.gz -b control_8_broden.bed -A > boundary_proc.tmp.bed

# bedtools subtract -a boundary_proc.tmp.bed -b exon.bed.gz -A -f 0.9 | \
bedtools subtract -a boundary.bed.gz -b control_8_broden.bed -A | \
    bedtools merge -i - | \
    awk '$1 ~ /^[0-9XY]+$/ {OFS="\t"; print $1, $2-10, $3}' > \
    boundary_filtered.hg19.grc.bed

rm -rf boundary.bed.gz*
rm -rf exon.bed.gz*
rm -rf control_8_broden.bed
rm -rf boundary_proc.tmp.bed

