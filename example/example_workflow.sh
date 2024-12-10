#!/bin/sh

# This file serves as an example workflow illustrating how and when to use msyd.
# It is a part of the msyd CI, and should pass so long as your system 


## Download three publicly available, high quality A. thaliana genomes

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

### Prepare them for running msyd

## rename them to shorter names
mv GCA_001651475.1_Ler_Assembly_genomic.fna ler.fna
mv GCA_028009825.1_Col-CC_genomic.fna col.fna
mv GCA_902460295.1_Arabidopsis_thaliana_Sha_genomic.fna sha.fna
mv GCA_024498555.1_ASM2449855v1_genomic.fna swe.fna

## filter out small scaffolds
grep -n -P ">" ./*.fna
# col and swe do not require truncating,
# for ler the small scaffolds start at line 1442097
head -n 1442096 ler.fna > ler.filtered.fna
mv ler.filtered.fna ler.fna
# and for sha at line 1480077
head -n 1480076 sha.fna > sha.filtered.fna
mv sha.filtered.fna sha.fna


### Generate inputs for msyd

## generate alignments to col-CC
mv col.fna ref.fna
minimap2 -cx asm10 --eqx ref.fna ler.fna > ler.paf
minimap2 -cx asm10 --eqx ref.fna sha.fna > sha.paf
minimap2 -cx asm10 --eqx ref.fna swe.fna > swe.paf

## run syri on the alignments
# make sure to pass --cigar and specify appropriate prefixes, so the msyd output is more easily interpretable
syri --nc 5 -F P --cigar --prefix ler -c ler.paf -r ref.fna -q ler.fna --lf ler.syri.log --samplename ler
syri --nc 5 -F P --cigar --prefix sha -c sha.paf -r ref.fna -q sha.fna --lf sha.syri.log --samplename sha
syri --nc 5 -F P --cigar --prefix swe -c swe.paf -r ref.fna -q swe.fna --lf swe.syri.log --samplename swe

## construct genomes.tsv file
# as msyd needs many input files, the paths are stored in a samplesheet
echo "#name\taln\tsyri\tvcf\tseq" > genomes.tsv
for f in *syri.out
do
	bs=$(basename $f syri.out)
	echo "$bs\t$bs.paf\t${bs}syri.out\t${bs}syri.vcf\t${bs}.fna" >> genomes.tsv
done

### run msyd to call multisynteny
msyd call -c 5 -i genomes.tsv -o athalianas.psf -m athalianas.vcf -r ref.fna

### work with the output

## export multisynteny on Chr3 for use in visualization/other software

# CP116283 is the id corresponding to chromosome 3 in Col-CC
# filter for multisynteny on this chromosome
msyd view -e "on CP116283.1" -i athalianas.psf -o athalianas-chr3.psf

# export to VCF; this could also be done in the command above
msyd view -i athalianas-chr3.psf -o athalianas-chr3-syn.vcf

## download 1001 genome project VCF, filter for small variants structurally conserved regions

curl https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-56/plants/variation/vcf/arabidopsis_thaliana/arabidopsis_thaliana.vcf.gz -o ensembl_athaliana.vcf.gz
gunzip ensembl_athaliana.vcf.gz

# change from ids to chr numbers, to match vcf nomenclature
sed -e s/CP116280\.1/Chr1/ -e s/CP116281\.1/Chr2/ -e s/CP116282\.1/Chr3/ -e s/CP116283\.1/Chr4/ -e s/CP116284\.1/Chr5/ athalianas.psf > athalianas-chrnames.psf

# filter for variants in coresyntenic regions!
msyd view -i athalianas-chrnames.psf -e "deg >= 3" -o coresynt-snvs.vcf --intersect ensembl_athaliana.vcf
