# Pansyn Filtering Language Reference

Pansyn uses user-defined logical expressions to allow filtering pansyntenic regions in a way that is more flexible and powerful than using simple `grep` commands.
The filtering is performed by evaluating the expression specified for each record individually, and adding that record to the output if the expression evaluates to `TRUE`.
The expression is a simple logical formula that can either be one of a number of properties of a record, or a combination of several such other expressions that can be nested to an arbitrary depth.
E.g. in order to filter for records on Chromosome 3 and containing (that is, being present in the genome of) the sample ler, the expression `(on Chr3
