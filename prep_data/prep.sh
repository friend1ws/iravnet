#! /bin/bash
#! -S /bin/bash
#$ -cwd
#$ -e log/ -o log/

###
# hg19 grc
annot_utils boundary boundary.bed.gz --donor_size 3,6 --acceptor_size 6,1 --grc

annot_utils exon --grc exon.bed.gz

zcat db/TCGA_742_4.bed.gz | awk '$1 ~ /^[0-9XY]+$/ {OFS="\t"; print $1, $2-10, $2+10}' > control_4_broden.bed

bedtools subtract -a boundary.bed.gz -b control_4_broden.bed -A > boundary_proc.tmp.bed

bedtools subtract -a boundary_proc.tmp.bed -b exon.bed.gz -A -f 0.9 | \
    bedtools merge -i - | \
    awk '$1 ~ /^[0-9XY]+$/ {OFS="\t"; print $1, $2-10, $3}' > \
    ../iravnet/data/boundary_proc.hg19.grc.bed

rm -rf boundary.bed.gz*
rm -rf exon.bed.gz*
rm -rf control_4_broden.bed
rm -rf boundary_proc.tmp.bed
###

###
# hg19
annot_utils boundary boundary.bed.gz --donor_size 3,6 --acceptor_size 6,1 

annot_utils exon exon.bed.gz

zcat db/TCGA_742_4.bed.gz | awk '$1 ~ /^[0-9XY]+$/ {OFS="\t"; print "chr"$1, $2-10, $2+10}' > control_4_broden.bed

bedtools subtract -a boundary.bed.gz -b control_4_broden.bed -A > boundary_proc.tmp.bed

bedtools subtract -a boundary_proc.tmp.bed -b exon.bed.gz -A -f 0.9 | \
    bedtools merge -i - | \
    awk '$1 ~ /^chr[0-9XY]+$/ {OFS="\t"; print $1, $2-10, $3}' > \
    ../iravnet/data/boundary_proc.hg19.bed

rm -rf boundary.bed.gz*
rm -rf exon.bed.gz*
rm -rf control_4_broden.bed
rm -rf boundary_proc.tmp.bed
###

###
# hg38 grc
annot_utils boundary boundary.bed.gz --donor_size 3,6 --acceptor_size 6,1 --genome_id hg38 --grc

annot_utils exon --genome_id hg38 --grc exon.bed.gz

zcat boundary.bed.gz > boundary_proc.tmp.bed

bedtools subtract -a boundary_proc.tmp.bed -b exon.bed.gz -A -f 0.9 | \
    bedtools merge -i - | \
    awk '$1 ~ /^[0-9XY]+$/ {OFS="\t"; print $1, $2-10, $3}' > \
    ../iravnet/data/boundary_proc.hg38.grc.bed

rm -rf boundary.bed.gz*
rm -rf exon.bed.gz*
rm -rf boundary_proc.tmp.bed
###

###
# hg19
annot_utils boundary boundary.bed.gz --donor_size 3,6 --acceptor_size 6,1 --genome_id hg38 

annot_utils exon --genome_id hg38 exon.bed.gz

zcat boundary.bed.gz > boundary_proc.tmp.bed

bedtools subtract -a boundary_proc.tmp.bed -b exon.bed.gz -A -f 0.9 | \
    bedtools merge -i - | \
    awk '$1 ~ /^chr[0-9XY]+$/ {OFS="\t"; print $1, $2-10, $3}' > \
    ../iravnet/data/boundary_proc.hg38.bed

rm -rf boundary.bed.gz*
rm -rf exon.bed.gz*
rm -rf boundary_proc.tmp.bed
###


