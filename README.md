# pansyn

## Installation
pansyn can be installed by a call to `setup.py` as `[sudo] python ./setup.py install [--user]`.
Requirements for running pansyn are `Cython`, `pysam`, `pandas` and `scipy` for the plotsr ordering functionality.

## Usage

Pansyn can be used as either a Python library or a CLI tool.
Pansyn can identify pansyntenic regions in a set of genomes aligned to a common reference, manipulate pansynteny callset files and determine a suitable ordering for plotting synteny between the genomes with `plotsr`.
The different functionalities of pansyn are separated into modules for the python library and subcommands for the CLI interface.

### `pansyn call`

This subcommand is used to call pansyntenic regions in a set of genomes.
It uses the Synteny Intersection Algorithm to identify both core and reference cross synteny, but can optionally be configured to call only core synteny (making it faster particularly when run with many genomes).
`pansyn call` expects a tab-separated input file (specified with `-i`) containing for each genome the name that should be used and the location of the alignment and SyRI.out file.
Both files should be using the same reference in all specified genomes.
An example input file might look something like the following:

```
\#name	aln	syri
an1	col_an1.bam	col_an1syri.out
c24	col_c24.bam	col_c24syri.out
eri	col_eri.bam	col_erisyri.out
```

`pansyn call` takes a number of optional CLI flags, described by calling `pansyn call -h`.

### `pansyn view`

The `view` subparser is used for operating on PFF files.
This command can also be used to convert output to VCF format (losing alignment information in the process); this requires a file containing the reference genomes.
`pansyn view` uses a custom expression language to allow complex filtering operations on PFF files.
Filtering expressions are constructed from two main building blocks: primitive and derived terms.
Primitive terms include e.g. filtering for the degree of a pansyn region (`deg >= 5`), for whether it exists in a query sequence (`contains sample01`) or its position on the reference (`in Chr1:x:200-500`).
There are two primitive terms that take multiple arguments (`contains all` and `contains any`).
In this case, the arguments are separated by a `,`. Any whitespace near the `,` is ignored.
Example: `contains any sample01, sample02,sample05`.
Derived terms alter or combine multiple other terms (these can in turn be primitive or derived).
If a derived term consists of multiple terms, these need to be put in brackets.
Example derived terms include `(deg >= 5) and (contains sample02)` or `not on Chr2`.

Internally, the expression is first parsed recursively into a lambda expression by `util.compile_filter`.
`util.filter_pansyns` then evaluates this for every pansyn region object, keeping only those for which the expression evaluates to `True`.
For users of the python API, `filter_pansyns` can also be called with any custom lambda expression.

For a complete reference of the filtering language, see `view_reference.md`

### `pansyn order`

`pansyn` can try to determine a suitable ordering for plotting a set of genomes once pansynteny has been called in that set.
The ordering functionality takes a pansynteny callset in PFF format as input (specified with `-i`) and prints the optimal ordering to stdout.
The organisms are given with the name that was specified in the input file to `pansyn call`.

The ordering is determined by performing single-linkage hierarchical clustering followed by leaf-order optimization (as implemented in `scipy`) according to the amount of shared pairwise pansynteny.
It is intended for within-species order determination; determining the optimal ordering for a particular region of the genome can be accomplished by filteringt he PFF file using `pansyn view`.
The algorithm may not produce suitable output if used on regions with little synteny.

## Example pansyn workflow

![Diagram illustrating an example workflow for using pansyn](https://github.com/schneebergerlab/pansyri/blob/leon/workflow.svg)

The first step in using `pansyn` is to assemble a number of query genomes and choose a reference.
Assemblies should be of high quality and be chromosome-scale; they may need to be scaffolded (using e.g. RagTag) before further processing.

In the next step, whole-genome alignments of the query sequences to the reference need to be computed.
`pansyn` can work with alignment files in SAM, BAM and PAF format.
Using minimap2, this step could be done using a script similar to the following:
```
#!/bin/bash

for seq in assemblies/*.fna; do
	bs=$(basename -s .fna $seq)
	minimap2 -cx asm5 --eqx ref.fna.gz $seq > alns/$bs.paf
done
```

After the alignments have been made, `syri` needs to be called for each of the queries:
```
#!/bin/bash

for seq in assemblies/*.fna; do
	bs=$(basename -s .fna $seq)
	syri --nc 5 -F P --cigar --dir syri --prefix $bs -c alns/$bs.paf -r ref.fna.gz -q $seq --lf $bs.syri.log
```

In preparation for running `pansyn call`, the input tsv needs to be generated.
An example TSV that could work with the code snippets above could look something like this:

```
#name	aln	syri
sample12	alns/sample12.paf	syri/sample12syri.out
sample20	alns/sample20.paf	syri/sample20syri.out
control3	alns/control3.paf	syri/control3syri.out
``` 

It could be generated manually, or with a bash script:
```
#!/bin/sh

echo "#name aln 	syri" > all.tsv
for f in syri/*syri.out
do
	bs=$(basename $f syri.out)
	echo "$bs	alns/$bs.paf	syri/${bs}syri.out" >> all.tsv
done
```

Finally, `pansyn` can be run on the dataset:
```
pansyn call -i all.tsv -o all.pff
```

Further processing can be done with `pansyn view`:


Useful scripts for using `pansyn` may be seen in the `hgsv\_scripts`, `drosophila\_scripts`, `ampsynth\_scripts` directories.
