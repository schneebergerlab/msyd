# The Population Syntenty File Format (PSF)
Updated to v0.3 file format 2024-09-11 Leon Rauschning
Updated to v0.2 file format 2024-05-22 Leon Rauschning
Created 2022-05-11 Leon Rauschning

## Premise

Finding synteny in multiple genomes at once presents a new way of working with the genomes of a population. To represent this novel kind of information, a new file format is needed; at the same time, as much compatibility as possible should be maintained with tools built for pairwise genome comparisons.
For this reason, in designing the Population Synteny File format we have chosen to retain some of the basic structure of the BED format and stick to a simple yet flexible tab-separated format.

This rough specification uses some terms specific to msyd; these are explained in glossary.md.

The information presented in PSF is richer than what msyd exports to VCF, in that it represents base-level alignments of multisyntenic regions and canrepresent synteny not found on a single linear reference.
To provide full compatibility with BED format, optionally even these annotations can be assigned a coordinate on a linear reference if specified; the region is then canonically placed at the earliest space it can be given the constraint of neighbouring regions.


## Spec
- A PSF file consists of a body and a header, and is tab-separated. A . represents a missing field, if it is allowed.
- The header starts with '\#' and contains the names of all organisms in the callset, separated by tabs
- The first three columns contain the genomic coordinates on a linear reference; if a sequence does not have a position on the reference, these may be annotated as missing values or as a single base at the earliest position allowed by the contraints of other regions present on the reference. The latter option is provided to maintain compatibility with BED and facilitate tabix indexing.
- The fourth column contains the ID of the multisyntenic annotation, labeling it as either coresyntenic (syntenic across all samples), merasyntenic (syntenic across only a subset of sequences) or private (not syntenic in any other sample). IDs present a way to consistently order all merasyntenic regions, even when they do not have a position on a linear reference.
- The fifth column contains the name of the sample chosen as representative for this multisyntenic region; it must be either present in the header annotation, or be 'ref' (the linear reference genome, by convention). The representative sample is the reference to which the alignments of the other samples are reported.
- The sixth to eighth columns store the position of the multisyntenic region on the representative sample.
- The following columns contain the position of each syntenic region in each sample and the CIGAR string of the alignment, separated by `,`. The CIGAR string may be elided for smaller file sizes and easier downstream processing, though further multisynteny identification may not work.
- The position is given as a range, consisting of an optional sample name specifier, a chromosome identifier and haplotype character separated by : from each other and the start and end positions separated by -. Both start and end positions are inclusive.
- A range is inverted if the start position is before the end position. In this case, any alignments referring to this range have this sequence reversed
- A range can contain whitespace characters that are not TAB at arbitrary positions; they are removed during parsing

## Example

CIGAR strings are shortened using [...]

\#CHR|START|END|ANN|REF|CHR|START|END|c24|eri|ler|sha
---
Chr1|124|530|MERASYN1|ref|.|.|.|.|Chr1:513-919,53=1X11=1X29=1X33=1X16=1X5=1X34=1X190=1X5=1X22=|.|.
Chr1|531|1084|MERASYN2|ref|.|.|.|.|Chr1:920-1458,44=1X27=1X94=1X12=1X43=1X136=1X133=15D26=1X4=1X2=1X3=1X5=|.|Chr1:1581-2126,44=1X21=1X5=1X22=1X71=1X12=1X43=1X136=1X173=8D1=1X8=1X
Chr1|1090|17001|CORESYN1|ref|.|.|.|Chr1:13-15966,1464=1X[...]|Chr1:1464-17369,1=2X[...]|Chr1:6-15923,5235=1I2344=[...]|Chr1:2132-18108,3=5I50=1X1=4D[...]
Chr1|17002|18729|MERASYN3|ref|.|.|.|Chr1:15967-17694,1X8=1X14=1X4=1X17=1X256=1X95=1X36=1X125=1I265=1X42=1X519=1X3=1X48=1X104=1D178=|.|Chr1:15924-17651,1055=1X672=|.
Chr1|18730|68988|CORESYN2|ref|.|.|.|Chr1:17695-67953,1793=10D3=1X[...]|Chr1:17378-67673,26=1X33=1X129=2I[...]|Chr1:17652-67896,229=1D1550=4D11696=[...]|Chr1:18117-67579,26=1X33=1X129=2I[...]
