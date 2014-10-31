repeatnet
=========

RepeatNet: an ab initio centromeric sequence detection algorithm

Additional tools that should be installed

* PHRAP for the assembly step
* TRF for calculating the consensus
* Graphviz for visualization   [ OPTIONAL ]


Compilation:
===========
    gcc repeatnet.c -o repeatnet -O3 -lm

input: FASTA format only.
The pairs of reads (preferably fosmid; if not plasmid) should be named
as:

>pair1.FORWARD.1
>pair1.REVERSE.1

>pair2.FORWARD.1
>pair2.REVERSE.1


the "pair1/pair2" part doesn't matter, but FORWARD/REVERSE part is used
to mark the pairs; of course the "pair1/pair2" part of the pairs must
match.

then run :
repeatnet -i test

this will generate a bunch of files:
test.h11.dump, test.h11.names, test.h11.winlog

then

repeatnet -loadwin test.h11.winlog -m -a 100 -c 100 -compare

-a 100 and -c 100  removes the vertices with less than 100 occurances in
the graph, the edges with <100 weight value (co-occurance frequency
between pairs of vertics).  Lower values will give more connected and
noisy graphs, higher values will "clean up" the graph more.

this will generate a *matrix and *viz file.  To visualize the graph, use
the graphviz package:

neato -Tps -o test.h11.winlog.h11.cut100.e10000.merged.eps test.h11.winlog.h11.cut100.e10000.merged.viz


then to divide into connected components:

repeatnet -loadmatrix test.h11.winlog.h11.cut100.e10000.merged.matrix  -loadnames test.h11.names  -components


pick the largest component first (by file size; or any other interesting-looking one
from the graph). in my test, it is component-30. then "encode" the kmers
back into numbers:

for i in `cut -f 1 component-30.txt`; do repeatnet -encode $i; done >
component-30.ids


this will calculate the hash values (or vertex id's) for the kmers.

then:

repeatnet -loadwin test.h11.winlog -loadnames test.h11.names -clones component-30.id


will generate a file component-30.ids.clones  that will have names of
the sequences (pairs) that are likely to contain satellite. This file
might be redundant, so run a "sort -u " on it.  Fetch the sequences back
from your original fasta file, and run phrap first;  then run TRF on the
contigs.





