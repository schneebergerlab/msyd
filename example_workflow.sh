#!/bin/sh

## Download some genomes

# download the Col-CC assembly
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCA_028009825.1/download?include_annotation_type=GENOME_FASTA&filename=GCA_028009825.1.zip" -H "Accept: application/zip"
# download the CRG Ler assembly
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCA_001651475.1/download?include_annotation_type=GENOME_FASTA&filename=GCA_001651475.1.zip" -H "Accept: application/zip"
# download the MPIPZ Sha assembly
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCA_902460295.1/download?include_annotation_type=GENOME_FASTA&filename=GCA_902460295.1.zip" -H "Accept: application/zip"
# download the GMI swedish assembly
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCA_024498555.1/download?include_annotation_type=GENOME_FASTA&filename=GCA_024498555.1.zip" -H "Accept: application/zip"

# unzip the datasets
unzip "./*.zip"
mv ncbi_dataset/data/*/*.fna ./

## Prepare them for running pasy

# rename them to shorter names
mv GCA_001651475.1_Ler_Assembly_genomic.fna ler.fna
mv GCA_028009825.1_Col-CC_genomic.fna col.fna
mv GCA_902460295.1_Arabidopsis_thaliana_Sha_genomic.fna sha.fna
mv GCA_024498555.1_ASM2449855v1_genomic.fna swe.fna

grep -n -P ">" ./*.fna
# col and swe do not require truncating,
# for ler the small scaffolds start at line 1442097
head -n 1442096 ler.fna > ler.filtered.fna
mv ler.filtered.fna ler.fna
# and for sha at line 1480077
head -n 1480076 sha.fna > sha.filtered.fna
mv sha.filtered.fna sha.fna


## Generate inputs for pasy

# generate alignments to col-CC
mv col.fna ref.fna
minimap2 -cx asm5 --eqx ref.fna ler.fna > ler.paf
minimap2 -cx asm5 --eqx ref.fna sha.fna > sha.paf
minimap2 -cx asm5 --eqx ref.fna swe.fna > swe.paf

# run syri on the alignments
syri --nc 5 -F P --cigar --prefix ler -c ler.paf -r ref.fna -q ler.fna --lf ler.syri.log --samplename ler
syri --nc 5 -F P --cigar --prefix sha -c sha.paf -r ref.fna -q sha.fna --lf sha.syri.log --samplename sha
syri --nc 5 -F P --cigar --prefix swe -c swe.paf -r ref.fna -q swe.fna --lf swe.syri.log --samplename swe

## construct genomes.tsv file
echo "#name\taln\tsyri\tvcf" > genomes.tsv
for f in *syri.out
do
	bs=$(basename $f syri.out)
	echo "$bs\t$bs.paf\t${bs}syri.out\t${bs}syri.vcf" >> genomes.tsv
done

# run pasy to call pansynteny
pasy call -i genomes.tsv -o athalianas.pff -m athalianas.vcf -r ref.fna

## work with the output

# CP116282 is the id corresponding to chromosome 3 in Col-CC
pasy view -e "on CP116283" -i athalianas.pff -o athalianas-chr3.pff

# convert to VCF for use in visualization/other software
pasy view -i athalianas-chr3.pff -o athalianas-chr3-syn.vcf

## download 1001 genome project VCF, filter for vars in pansyntenic regions

curl https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-56/plants/variation/vcf/arabidopsis_thaliana/arabidopsis_thaliana.vcf.gz -o ensembl_athaliana.vcf.gz
gunzip ensembl_athaliana.vcf.gz

# change from ids to chr numbers, to match vcf nomenclature
sed -e s/CP116280.1/1/ -e s/CP116281.1/2/ -e s/CP116282.1/3/ -e s/CP116283.1/4/ -e s/CP116284.1/5/ athalianas.pff > athalianas-chrnames.pff

# filter for variants in pansyntenic regions!
pasy view -i athalianas-chrnames.pff -e "deg >= 3" -o pansynt-vars.vcf --intersect ensembl_athaliana.vcf
