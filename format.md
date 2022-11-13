# The PanSyRI File Format

## Basic premise
Multigenomic structure is more complex and less straightforward to fit into categories than genomic structure across only two genomes.
This necessitates a more flexible file format, able to account for complex structural landscapes.
This format can store either just the annotation produced by PanSyRI or the annotation and all or part of the corresponding alignments.
It is intended to be straightforward both for human and automatic parsing.

## Spec
- A PFF file consists of a body and a header, and is tab-separated
- the header starts with '\#ANN' and contains the names of all organisms in the callset, separated by tabs
- If a reference was used, it is considered part of the callset and given first in the header. It should be named 'ref'.
- Each row stores one annotation across the set of organisms
- An annotation consists of an annotation type and a list of ranges for each organism, separated by ;
- A range consists of a chromosome number and haplotype character separated by : from each other and the start and end positions separated by -, and optionally a list of alignment specifiers, separated by ,s
- x represents an unknown haplotype; X represents a locus in all haplotypes
- A range can contain whitespace characters that are not TAB at arbitrary positions; they are removed during parsing
- An alignment is specified by first providing the organism, then the 1-based index of the range in that organism's field to which this range has been aligned and finally the alignment in CIGAR format, all separated by :
- If a reference is specified and the alignment is to the reference and the reference field only has one range, just the CIGAR string can be given instead.
- For any given CIGAR string, the reference is the organism corresponding specified by the alignment specifier and the query is the column the CIGAR string occurs at.
- A range is inverted if the start position is lower than the end position.
- In this case, any alignments referring to this range have this sequence reversed
– the SNP annotation is treated in a special manner:
	- A SNP is stored as a position, that is a range with no end
	- SNP annotations cannot store CIGAR strings; after a , there may only be one character storing the genotype at that locus


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

## Probable current solutions
https://github.com/lh3/gfatools/blob/master/doc/rGFA.md
