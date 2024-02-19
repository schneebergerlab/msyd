# An example msyd Workflow

This file describes the steps for an example workflow using msyd, described in the image below.
See also `example_workflow.sh` for a condensed shell script of the same workflow that can be run on the command line.


![Diagram illustrating an example workflow for using msyd](https://github.com/schneebergerlab/msyd/blob/master/img/workflow.svg)

## Running msyd

In order to use `msyd`, high-quality chromosome-scale genomes are required.
Scaffolding using e.g. RagTag may be necessary.
For this example workflow, we will be downloading some A. thaliana assemblies from the GenBank database:

```
# download the Col-CC assembly
$ curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCA_028009825.1/download?include_annotation_type=GENOME_FASTA&filename=GCA_028009825.1.zip" -H "Accept: application/zip"
# download the CRG Ler assembly
$ curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCA_001651475.1/download?include_annotation_type=GENOME_FASTA&filename=GCA_001651475.1.zip" -H "Accept: application/zip"
# download the MPIPZ Sha assembly
$ curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCA_902460295.1/download?include_annotation_type=GENOME_FASTA&filename=GCA_902460295.1.zip" -H "Accept: application/zip"
# download the GMI swedish assembly
$ curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCA_024498555.1/download?include_annotation_type=GENOME_FASTA&filename=GCA_024498555.1.zip" -H "Accept: application/zip"

# unzip the datasets
$ unzip "./*.zip"
$ mv ncbi_dataset/data/*/*.fna ./

# rename them to shorter names
$ mv GCA_001651475.1_Ler_Assembly_genomic.fna ler.fna
$ mv GCA_028009825.1_Col-CC_genomic.fna col.fna
$ mv GCA_902460295.1_Arabidopsis_thaliana_Sha_genomic.fna sha.fna
$ mv GCA_024498555.1_ASM2449855v1_genomic.fna swe.fna
```

Some of these assemblies still contain short scaffolds not corresponding to any of the chromosomes:

```
$ grep -P ">" ./*.fna
ler.fna:>CM004359.1 Arabidopsis thaliana ecotype Landsberg erecta chromosome 1, whole genome shotgun sequence
ler.fna:>CM004360.1 Arabidopsis thaliana ecotype Landsberg erecta chromosome 2, whole genome shotgun sequence
ler.fna:>CM004361.1 Arabidopsis thaliana ecotype Landsberg erecta chromosome 3, whole genome shotgun sequence
ler.fna:>CM004362.1 Arabidopsis thaliana ecotype Landsberg erecta chromosome 4, whole genome shotgun sequence
ler.fna:>CM004363.1 Arabidopsis thaliana ecotype Landsberg erecta chromosome 5, whole genome shotgun sequence
ler.fna:>LUHQ01000006.1 Arabidopsis thaliana scaffold15_Contig142, whole genome shotgun sequence
ler.fna:>LUHQ01000007.1 Arabidopsis thaliana scaffold15_Contig624, whole genome shotgun sequence
ler.fna:>LUHQ01000008.1 Arabidopsis thaliana scaffold18_size294915, whole genome shotgun sequence
ler.fna:>LUHQ01000009.1 Arabidopsis thaliana scaffold24_size307384, whole genome shotgun sequence
ler.fna:>LUHQ01000010.1 Arabidopsis thaliana scaffold26_size238942, whole genome shotgun sequence
ler.fna:>LUHQ01000011.1 Arabidopsis thaliana scaffold27_size282142, whole genome shotgun sequence
ler.fna:>LUHQ01000012.1 Arabidopsis thaliana scaffold29_size187832, whole genome shotgun sequence
[snip]
```

As the scaffolds are ordered by size, the chromosomes will always be first.
We can use `grep` again with the `-n` option to find the starting line of the first non-chromosomal scaffold and then use the `head` command to filter them from the FASTA:

```
$ grep -n -P ">" ./*.fna
# col and swe do not require truncating,
# for ler the small scaffolds start at line 1442097
# and for sha at line 1480077
$ head -n 1442096 ler.fna > ler.filtered.fna
$ mv ler.filtered.fna ler.fna
$ head -n 1480076 sha.fna > sha.filtered.fna
$ mv sha.filtered.fna sha.fna
```


In the next step, whole-genome alignments of the query sequences to the reference need to be computed.
`msyd` can work with alignment files in SAM, BAM and PAF format.
We choose Columbia as a reference, and align the sequences using `minimap2`:

```
$ mv col.fna ref.fna
$ minimap2 -cx asm5 --eqx ref.fna ler.fna > ler.paf
$ minimap2 -cx asm5 --eqx ref.fna sha.fna > sha.paf
$ minimap2 -cx asm5 --eqx ref.fna swe.fna > swe.paf
```

After the alignments have been made, `syri` needs to be run on each of the alignments:

```
$ syri --nc 5 -F P --cigar --prefix ler -c ler.paf -r ref.fna -q ler.fna --lf ler.syri.log --samplename ler
$ syri --nc 5 -F P --cigar --prefix sha -c sha.paf -r ref.fna -q sha.fna --lf sha.syri.log --samplename sha
$ syri --nc 5 -F P --cigar --prefix swe -c swe.paf -r ref.fna -q swe.fna --lf swe.syri.log --samplename swe
```

In preparation for running `msyd call`, the input tsv needs to be generated.
An example TSV for the callset prepared in this workflow could look like this:

```
$ cat genomes.tsv
#name   aln syri    vcf seq
ler	ler.paf	lersyri.out    lersyri.vcf  ler.fna
sha	sha.paf	shasyri.out    shasyri.vcf  sha.fna
swe	swe.paf	swesyri.out    swesyri.vcf  swe.fna
``` 

It could be typed out manually, or automatically generated with a bash script:

```
#!/bin/sh
# \t is the tab character, escaped here because it is converted to spaces by some markdown renderers
echo "#name\taln\tsyri\tvcf" > genomes.tsv
for f in *syri.out
do
	bs=$(basename $f syri.out)
	echo "$bs\t$bs.paf\t${bs}syri.out\t${bs}syri.vcf\t${bs}.fna" >> genomes.tsv
done
```

After generating the requisite input, `msyd call` can be run to generate the pansynteny callset:

```
$ msyd call -i genomes.tsv -o athalianas.pff -m athalianas.vcf -r ref.fna
```

## Working with the msyd Output

If we are only interested in pansynteny on Chromosome 3, we can filter the .PFF file calling msyd view.
In this case, the filtering could also be done using `grep` or `awk`, but `msyd view` is advantageous for more complex filtering.
A simple expression-based language is used for complex filtering instructions in `msyd`.
It is described in `view_reference.md`.

```
# CP116282 is the id corresponding to chromosome 3 in Col-CC
$ msyd view -e "on CP116283.1" -i athalianas.pff -o athalianas-chr3.pff
```

If we want to have the pansynteny annotations in VCF format, we can just change the file extension passed to `msyd view` or pass the `--ovcf` flag.

```
$ msyd view -i athalianas-chr3.pff -o athalianas-chr3-syn.vcf
```

We can also use our pansynteny callset to select variants in pansyntenic regions from a population VCF for which we do not have a assembled genomes.
As an example, let's download the 1001 genomes project VCF from the EBI:

```
$ curl https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-56/plants/variation/vcf/arabidopsis_thaliana/arabidopsis_thaliana.vcf.gz -o ensembl_athaliana.vcf.gz
$ gunzip ensembl_athaliana.vcf.gz
```

Before we can filter the file, we need to change the contig names used in `athalianas.pff` to chromosome names (at least for the reference).
We can do this with using `sed` to replace contig IDs with the corresponding chromosome number:

```
$ sed -e s/CP116280.1/1/ -e s/CP116281.1/2/ -e s/CP116282.1/3/ -e s/CP116283.1/4/ -e s/CP116284.1/5/ athalianas.pff > athalianas-chrnames.pff
```


For this, we can use the `--intersect` flag.
We can at the same time specify a filtering option using the `-e` flag.
Filtering will be performed before the VCF subsetting – this can be used to select only core-syntenic regions:

```
$ msyd view -i athalianas-chrnames.pff -e "deg >= 3" -o pansynt-vars.vcf --intersect ensembl_athaliana.vcf
```
