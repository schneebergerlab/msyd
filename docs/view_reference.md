# msyd Filtering Language Reference

msyd uses user-defined logical expressions to allow filtering pansyntenic regions in a way that is more flexible and powerful than using simple `grep` commands.
The filtering is performed by evaluating the expression specified for each record individually, and adding that record to the output if the expression evaluates to `TRUE`.
The expression is a simple logical formula that can either be one of a number of properties of a record, or a combination of several such other expressions that can be nested to an arbitrary depth.
E.g. in order to filter for records on Chromosome 3 and containing (that is, being present in the genome of) the sample ler, the expression `(on Chr3) and (contains ler)` can be used.
Many expressions can also be abbreviated (for example the above expression can also be written as `(on Chr3) & (cont ler)`).
Keywords are not case sensitive, so `AND`, `And` or even `aNd` also works, but sample names need to match the case that was specified in the TSV file.
See below for a full table of the expressions and how they can be shortened.
Expressions can be nested to an arbitrary depth.
So if we wanted to modify the expression above to also allow records that contain sha, but not ler, we could write it as `(on Chr3) and ((cont ler) or (cont sha))`.

expression | what it does
-- | --
not X, !X | accepts a record if X evaluates to False, otherwise declines
(X) and, & (Y) | accepts a record if both X and Y evaluate to True
(X) or, \| (Y) | accepts a record if any of X and Y evaluate to True
(X) xor, ^ (Y) | accepts a record if one, but not both of X and Y evaluate to True
True | always accepts a record
False | never accepts a record
(X) | same as X (in case there are extra brackets lying around)
deg >=/<= INT | true if the record has a degree >= or <= than INT (which needs to be a number)
len >=/<= INT | true if the record has a length <= or >= than INT on the reference
contains, cont ORG | true if the record has an entry for the organism ORG. ORG needs to match case with the organism name specified in the TSV!
containsall, contall ORG1, ORG2, ... | true if the record has an entry for all of ORG1, ORG2, etc.
containsany, contany ORG1, ORG2, ... | true if the record has an entry for at least one of ORG1, ORG2, etc.
on CHR | true if the record lies on Chromosome CHR on the reference (CHR needs to match the case used in the PFF File)
in CHR:START-END | true if the records location on the reference lies on Chromosome CHR between START and END (both integers, CHR needs to match the case used in the PFF File).




