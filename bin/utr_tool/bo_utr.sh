#!/bin/bash
set -e

#Inputs:
# $1 - Genome sequence fasta file (genome.fa)
# $2 - Genome feature / annotation file (annot.gff3)
# $3 - Gene ID prefix to look for in annot.gff3
#      This is a string such as "gene_id=" or "ID=gene:" which prefixes the gene ID in column 9 of a gff3 file)
# $4 - Length of promoter sequence to take (upstream of TSS)
# $5 - Toggle whether to remove promoter overlapping bits with gene sequences
# $6 - Toggle whether to include 5' UTR sequence


# Prepare workspace
if [ ! -d "works" ]; then
	mkdir works
fi

# if [ ! -d "outputs" ]; then
# 	mkdir outputs
# fi

mv $1 works/genome.fa
mv $2 works/annot.gff3

cd works

# 

#start off by filtering the .gff3 to gene lines only
# cp $2 works/annot.gff3
grep -P '\tgene\t' annot.gff3 > genelines.gff3

#strip the potential FASTA line breaks. creates genome_stripped.fa
# cp $1 genome.fa
python3 ../scripts/strip_newlines.py # genome.fa -> genome_stripped.fa
rm genome.fa

#create the .genome file
samtools faidx genome_stripped.fa # creates .fai file
cut -f 1-2 genome_stripped.fa.fai > bedgenome.genome

#parse up the .bed for promoter extraction
#the python script takes the genelines.gff3 file and makes a genelines.bed out of it
python3 ../scripts/parse_genelines.py $3 # -> genelines.bed
python3 ../scripts/parse_utrs_bo.py

#create universe
cut -f 4 genelines.bed > universe.txt

bedtools flank -l $4 -r 0 -s -i genelines.bed -g bedgenome.genome > promoters.bed

#possibly add 5' UTR
if [ $6 == '--UseUTR' ]
	then
		python3 ../scripts/parse_utrs.py
fi

#remove overlapping promoter chunks
if [ $5 == '--NoOverlap' ]
	then
		#genelines2.bed is a CDS-only region file created in parse_utrs_bo.py
		python3 ../scripts/subtract.py --Input promoters.bed --Remove genelines2.bed --Output promoters2.bed --Kick
		python3 ../scripts/subtract.py --Input utr5.bed --Remove genelines2.bed --Output utr52.bed
		mv promoters2.bed promoters.bed
		mv utr52.bed utr5.bed
fi
#assess integrity
python3 ../scripts/assess_integrity.py

echo "UTR version 1.1.1"


bedtools getfasta -fi genome_stripped.fa -bed promoters.bed -s -fo promoters.fa -name
#this results in some really crappy nomenclature for gene names
#so let's make promoters.fa ourselves
#python3 ./scripts/parse_promoters.py

#there's a lot of intermediate files that need blanking. leave promoters.fa
# rm genelines.bed
# rm genelines.gff3
# rm genome.fa
# rm genome_stripped.fa
# rm genome_stripped.fa.fai
# rm promoters.bed

# rm universe.txt