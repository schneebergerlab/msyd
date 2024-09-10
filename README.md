# msyd

msyd is still under active development, so expect some bugs and changes!
If in doubt about the behaviour of msyd or how it might change, feel free to reach out by opening an issue!

## Changelog
2024-08-23: As of now, the realignment step requires a patched version of minimap2 that exposes some configuration options not exposed in the default version. Compare also lh3/minimap2#1240. The patched version is available [here](https://github.com/lrauschning/minimap2/tree/mappy) temporarily.

## Installation
msyd can be installed by a call to `setup.py` as `python ./setup.py install [--user]` or with `pip install .` in the base directory of the repository.
Requirements for running msyd are `python >= 3.8`, `Cython`, `pysam >= 0.21.0`, `pandas`, `numpy`, `scipy`, `mappy` (minimap2; see changelog above), `intervaltree` and `syri >= 1.6.5`, as well as a C++ compiler supporting OpenMP.
The requirements are listed in `requirements.txt` and can be installed via conda or another package manager from this file.

## Usage

msyd can be used as either a Python library or a CLI tool.
msyd can identify multisyntenic regions in a set of genomes aligned to a common reference, manipulate multisynteny callset files and determine a suitable ordering for plotting synteny between the genomes with `plotsr`.
The different functionalities of msyd are separated into modules for the python library and subcommands for the CLI interface.
Many convenient shorthand functions are defined in the `util` module.

### `msyd call`

This subcommand is used to call multisyntenic regions in a set of genomes.
It uses the Synteny Intersection Algorithm to identify both core and reference merasynteny, but can optionally be configured to call only coresynteny (making it faster particularly when run with many genomes).
`msyd call` expects a tab-separated input file (specified with `-i`) containing for each genome the name that should be used and the location of the alignment, SyRI.out and vcf file (make sure the SyRI.out file was generated without the `--maxsize` parameter, otherwise large INDELs will be missing!).
Optionally, the VCF might be left out, then `msyd` will automatically look for the syri-generated VCF file corresponding to the syri.out file.
Both files should be using the same reference in all specified genomes.
An example input file might look something like the following:

```
$cat genomes.tsv
#name	aln	syri	vcf	seq
an1	col_an1.bam	col_an1syri.out	col_an1.vcf	an1.fa.gz
c24	col_c24.bam	col_c24syri.out	col_c24.vcf	c24.fa.gz
eri	col_eri.bam	col_erisyri.out	col_eri.vcf	eri.fa.gz
```

Multisynteny can then be called with `msyd`:
```
$msyd call -i genomes.tsv -o output_threesamples.pff
```


`msyd call` takes a number of optional CLI flags, described by calling `msyd call -h`.
If the -m output is specified with an output vcf, msyd will merge VCF files of each organism againt the reference into one large, multi-genomic VCF File.
This functionality only works for VCF files generated with a recent (>=1.6.3) version of SyRI.
Also be sure to pass the `--samplename` parameter to syri, preferably using the same name as specified in the .tsv.
In the call above, the VCF merging functionality can then be enabled:
```
$msyd call -i genomes_with_vcf.tsv -o threesamples.pff -r col.fa.gz -m threesamples.vcf
```

`msyd call` generates output in Population synteny File Format (PFF).
PFF is a tab-separated format for storing multisynteny annotations and alignments.
A specification and small example can be seen in `format.md`.
`msyd view` can be used to export the synteny annotations in a PFF file into a VCF for use with other tools, but `msyd` cannot read multisynteny information from VCF files.

### `msyd view`

The `view` subparser is used for operating on PFF files.
This command can also be used to convert output to VCF format (losing alignment information in the process); this requires a file containing the reference genomes.
`msyd view` uses a custom expression language to allow complex filtering operations on PFF files.
Filtering expressions are constructed from two main building blocks: primitive and derived terms.
Primitive terms include e.g. filtering for the degree of a msyd region (`deg >= 5`), for whether it exists in a query sequence (`contains sample01`) or its position on the reference (`in Chr1:x:200-500`).
There are two primitive terms that take multiple arguments (`contains all` and `contains any`).
In this case, the arguments are separated by a `,`. Any whitespace near the `,` is ignored.
Example: `contains any sample01, sample02,sample05`.
Derived terms alter or combine multiple other terms (these can in turn be primitive or derived).
If a derived term consists of multiple terms, these need to be put in brackets.
Example derived terms include `(deg >= 5) and (contains sample02)` or `not on Chr2`.

Internally, the expression is first parsed recursively into a lambda expression by `util.compile_filter`.
`util.filter_multisyns` then evaluates this for every msyd region object, keeping only those for which the expression evaluates to `True`.
For users of the python API, `filter_multisyns` can also be called with any custom lambda expression.

For a complete reference of the filtering language, see `view_reference.md`

If the `--intersect` option is passed with a VCF file, a VCF file containing only variants annotated in multisyntenic regions will be written to the file specified in `-o`.
If `-e` is also specified, the filtering will be performed before the intersection.
An example usecase would be to only get multisyntenic SNPs from Chr. 8:
```
msyd view -e 'on Chr8' -i calls.pff -e --intersect snps.vcf
```

### `msyd order`

`msyd` can try to determine a suitable ordering for plotting a set of genomes once multisynteny has been called in that set.
The ordering functionality takes a multisynteny callset in PFF format as input (specified with `-i`) and prints the optimal ordering to stdout.
The organisms are given with the name that was specified in the input file to `msyd call`.

The ordering is determined by performing single-linkage hierarchical clustering followed by leaf-order optimization (as implemented in `scipy`) according to the amount of shared pairwise multisynteny.
It is intended for within-species order determination; determining the optimal ordering for a particular region of the genome can be accomplished by filtering the PFF file using `msyd view`.
The algorithm may not produce suitable output if used on regions with little synteny.

## Example msyd workflow

![Diagram illustrating an example workflow for using msyd](https://github.com/schneebergerlab/msyd/blob/master/img/workflow.svg)

Useful scripts for using `msyd` may be seen in the `hgsv\_scripts`, `drosophila\_scripts`, `ampsynth\_scripts` directories.

# Contact

Please report any crashes, weird output, suggestions or otherwise in a GitHub issue!
