# An example pasy Workflow
![Diagram illustrating an example workflow for using pasy](https://github.com/schneebergerlab/pasy/blob/leon/workflow.svg)

In order to use `pasy`, high-quality chromosome-scale genomes are required.
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

# unzip the datasets, move to seqs folder
$ unzip ./*.zip
$ mkdir seqs
$ mv ncbi_datasets/data/*/*.fna seqs/

# remove the zipped files
$ rm -r ncbi_datasets
$ rm *.zip

# rename them to shorter names
$ cd seqs
$ mv GCA_001651475.1_Ler_Assembly_genomic.fna ler.fna
$ mv GCA_028009825.1_Col-CC_genomic.fna col.fna
$ mv GCA_902460295.1_Arabidopsis_thaliana_Sha_genomic.fna sha.fna
$ mv GCA_024498555.1_ASM2449855v1_genomic.fna swe.fna
$ cd ..
```

Some of these assemblies still contain short scaffolds not corresponding to any of the chromosomes:

```
$ grep -P ">" seqs/*.fna
seqs/ler.fna:>CM004359.1 Arabidopsis thaliana ecotype Landsberg erecta chromosome 1, whole genome shotgun sequence
seqs/ler.fna:>CM004360.1 Arabidopsis thaliana ecotype Landsberg erecta chromosome 2, whole genome shotgun sequence
seqs/ler.fna:>CM004361.1 Arabidopsis thaliana ecotype Landsberg erecta chromosome 3, whole genome shotgun sequence
seqs/ler.fna:>CM004362.1 Arabidopsis thaliana ecotype Landsberg erecta chromosome 4, whole genome shotgun sequence
seqs/ler.fna:>CM004363.1 Arabidopsis thaliana ecotype Landsberg erecta chromosome 5, whole genome shotgun sequence
seqs/ler.fna:>LUHQ01000006.1 Arabidopsis thaliana scaffold15_Contig142, whole genome shotgun sequence
seqs/ler.fna:>LUHQ01000007.1 Arabidopsis thaliana scaffold15_Contig624, whole genome shotgun sequence
seqs/ler.fna:>LUHQ01000008.1 Arabidopsis thaliana scaffold18_size294915, whole genome shotgun sequence
seqs/ler.fna:>LUHQ01000009.1 Arabidopsis thaliana scaffold24_size307384, whole genome shotgun sequence
seqs/ler.fna:>LUHQ01000010.1 Arabidopsis thaliana scaffold26_size238942, whole genome shotgun sequence
seqs/ler.fna:>LUHQ01000011.1 Arabidopsis thaliana scaffold27_size282142, whole genome shotgun sequence
seqs/ler.fna:>LUHQ01000012.1 Arabidopsis thaliana scaffold29_size187832, whole genome shotgun sequence
[snip]
```

As the scaffolds are ordered by size, the chromosomes will always be first.
We can use `grep` again with the `-n` option to find the starting line of the first non-chromosomal scaffold and then use the `head` command to filter them from the FASTA:

```
$ grep -n -P ">" seqs/*.fna
# col and we do not require truncating,
# for ler the small scaffolds start at line 1442097
# and for sha at line 1480077
$ head -n 1442096 seqs/ler.fna > seqs/ler.filtered.fna
$ mv seqs/ler.filtered.fna seqs/ler.fna
$ head -n 1480076 seqs/sha.fna > seqs/sha.filtered.fna
$ mv seqs/sha.filtered.fna seqs/sha.fna
```


In the next step, whole-genome alignments of the query sequences to the reference need to be computed.
`pasy` can work with alignment files in SAM, BAM and PAF format.
We choose Columbia as a reference, and align the sequences using `minimap2`:

```
$ mv seqs/col.fna ./ref.fna
$ mkdir alns
$ minimap2 -cx asm5 --eqx ref.fna seqs/ler.fna > alns/ler.paf
$ minimap2 -cx asm5 --eqx ref.fna seqs/sha.fna > alns/sha.paf
$ minimap2 -cx asm5 --eqx ref.fna seqs/swe.fna > alns/swe.paf
```

After the alignments have been made, `syri` needs to be run on each of the alignments:

```
$ syri --nc 5 -F P --cigar --dir syri --prefix ler -c alns/ler.paf -r ref.fna -q seqs/ler.fna --lf ler.syri.log
$ syri --nc 5 -F P --cigar --dir syri --prefix sha -c alns/sha.paf -r ref.fna -q seqs/sha.fna --lf sha.syri.log
$ syri --nc 5 -F P --cigar --dir syri --prefix swe -c alns/swe.paf -r ref.fna -q seqs/swe.fna --lf swe.syri.log
```

In preparation for running `pasy call`, the input tsv needs to be generated.
An example TSV that could work with the code snippets above could look something like this:

```
$ cat genomes.tsv
#name	aln	        syri                vcf
ler	alns/ler.paf	syri/lersyri.out    syri/lersyri.vcf
sha	alns/sha.paf	syri/shasyri.out    syri/shasyri.vcf
swe	alns/swe.paf	syri/swesyri.out    syri/swesyri.vcf
``` 

It could be typed out manually, or automatically generated with a bash script:

```
#!/bin/sh

echo "#name aln syri    vcf" > all.tsv
for f in syri/*syri.out
do
	bs=$(basename $f syri.out)
	echo "$bs	alns/$bs.paf	syri/${bs}syri.out  syri/${bs}syri.vcf" >> all.tsv
done
```

After generating the requisite input, `pasy call` can be run to generate the pansynteny callset:

```
$ pasy call -i genomes.tsv -o athalianas.pff -m athalianas.vcf -r ref.fna
```

Further processing can be done with `pasy view`:


