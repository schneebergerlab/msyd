# File Format draft

## Basic Ideas
- Store variants as a list of ranges/positions for each organism
- A range consists of a chromosome number, haplotype character, start and end positions separated by a colon, and a ?CIGAR string? detailing the match to a reference, separated by a comma
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
SYN|	5X10:100,90=|	5X10:100,90=|	5a10:100,45=5X30=10D;5b10:100,90=
TRANS|	8X20:30,10=;9a440:450,4=2X4=|	8X20:30,10=|	8X20:30,8=2D
SNP|	12X400:400,1=|	12X400:400,1X|	12a400:400,1=;12b400:4001X
INV|	2X100:200,100=|	2X200:100,100=|	2X100:200
DEL|	4X10:20,10=|	|	4a10:20,10=
INS|	|	10b30:50,20I|	10X30:50,20I
DUP|	5X20:30,10=|	5X20:30,10=;5a30:40,10I|	5X20:30,10=
