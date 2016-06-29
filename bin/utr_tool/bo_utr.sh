#!/bin/bash
set -e

#inputs:
# $1 - genome
# $2 - gff3
# $3 - gff3 header to use (TAIR 'gene_id=')
# $4 - toggle whether to remove promoter overlapping bits with gene sequences
# $5 - toggle whether to include 5' UTR sequence

#start off by filtering the .gff3 to gene lines only
cp $2 annot.gff3
grep -P '\tgene\t' annot.gff3 > genelines.gff3

#strip the potential FASTA line breaks. creates genome_stripped.fa
cp $1 genome.fa
python3 /scripts/strip_newlines.py

#create the .genome file
samtools faidx genome_stripped.fa
cut -f 1-2 genome_stripped.fa.fai > bedgenome.genome

#parse up the .bed for promoter extraction
python3 /scripts/parse_genelines.py $3
#the python script takes the genelines.gff3 file and makes a genelines.bed out of it
bedtools flank -l $4 -r 0 -s -i genelines.bed -g bedgenome.genome > promoters.bed
#remove overlapping promoter chunks
if [ $4 == '--NoOverlap' ]
	then
		bedtools subtract -a promoters.bed -b genelines.bed > promoters2.bed
		mv promoters2.bed promoters.bed
fi
#possibly add 5' UTR
if [ $5 == '--UseUTR' ]
	then
		python3 /scripts/parse_utrs.py
fi
bedtools getfasta -fi genome_stripped.fa -bed promoters.bed -s -fo promoters_rough.fa
#this results in some really crappy nomenclature for gene names
#so let's make promoters.fa ourselves
python3 /scripts/parse_promoters.py

#there's a lot of intermediate files that need blanking. leave promoters.fa
rm genelines.bed
rm genelines.gff3
rm genome.fa
rm genome_stripped.fa
rm genome_stripped.fa.fai
rm promoters.bed
rm promoters.bg
rm promoters_rough.fa
rm promoter_lengths.txt