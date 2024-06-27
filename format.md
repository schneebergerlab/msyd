# The Population syntenty File Format (PFF)

## Premise
Finding synteny in multiple genomes at once presents novel challenges, as not all syntenic sequences may be present on a linear reference.
To usefully represent this information, a new file format is needed; at the same time, as much compatibility as possible should be maintained with tools built for pairwise genome comparisons.
For this reason, PFF retains some of the basic structure of the BED format.
The information presented in PFF is richer than what msyd exports to VCF, in particular it is able to contain synteny not present in the linear reference.
For compatibility with BED format, a position on a linear reference is reported for all annotations, even those not present on the reference; in these cases, they are placed at the first position they can be at given the constraint of the multisyntenic regions surrounding them.

## Spec
- A PFF file consists of a body and a header, and is tab-separated
- The header starts with '\#ANN' and contains the names of all organisms in the callset, separated by tabs
- The first three columns contain the genomic coordinates on a linear reference; if a sequence does not have a position on the reference, it is annotated as a single base somewhere between the two nearest regions that constrain its position on the reference.
- The fourth column contains the ID of each syntenic region, labeling it as either coresyntenic (syntenic across all samples) or merasyntenic (syntenic across only a subset of sequences). IDs present a way to consistently order all merasyntenic regions, even when they do not have a position on a linear reference; if possible, order of genomic coordinates on a reference are preserved.
- The following columns contain the position of each syntenic region in each sample, as well as the sample this region was aligned against to make the call and optionally the CIGAR string of the alignment, separated by `,`.
- The position is given as a range, consisting of an optional sample name specifier, a chromosome number and haplotype character separated by : from each other and the start and end positions separated by -
- x represents an unknown haplotype; X represents a locus in all haplotypes
- A range is inverted if the start position is before the end position. In this case, any alignments referring to this range have this sequence reversed
- A range can contain whitespace characters that are not TAB at arbitrary positions; they are removed during parsing
- Merasyntenic regions, in particular those not present in the reference may be collapsed into one record by separating the annotations in each column with a `;`

## Example

ANN|	Ref|	Qry1|	Qry2
-|	-|	-|	-
SYN|	5:X:10-100,90=|	5:X:10-100,90=|	5:a:10-100,45=5X30=10D;5:b:10-100,90=
TRANS|	8:X:20-30,10=;9:a:440-450,4=2:X:4=|	8:X:20-30,10=|	8:X:20-30,8=2D
SNP|	12:X:400-400,1=|	12:X:400-400,1:X:|	12:a:400-400,1=;12:b:400-4001:X:
INV|	2:X:100-200,100=|	2:X:200-100,100=|	2:X:100-200
DEL|	4:X:10-20,10=|	|	4:a:10-20,10=
INS|	|	10:b:30-50,20I|	10:X:30-50,20I
DUP|	5:X:20-30,10=|	5:X:20-30,10=;5:a:30-40,10I|	5:X:20-30,10=

