# File Format draft

## Basic Ideas
- Store variants as a list of ranges/positions for each organism
- A range consists of a chromosome number and haplotype character separated by : from each other and the start and end positions separated by -, and optionally a ?CIGAR string? detailing the match to a reference, separated by a comma from the end position
- x represents an unknown haplotype; X represents a position in all haplotypes
- A range is inverted if the start position is lower than the end position. In this case, the CIGAR string represents a reversed alignment
- A range can contain whitespace characters that are not TAB at arbitrary positions; they are removed during parsing
- Each row stores one region across the organisms
- Each row starts with a column containing an annotation, i.e. SYN, DEL, INS
- columns are separated by a TAB character
- columns contain a list of ranges separated by semicolons
- a SNP is stored as a position, that is a range with length one (or without end?)
- TODO how to define the reference for CIGAR strings?

## Example

Annotation|	Ref|	Qry1|	Qry2
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
