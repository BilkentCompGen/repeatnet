# TODO

* Hash table and the hashtest function seem to be for finding the read pairs if the input is not correctly interleaved. We can use this only if the input is not correctly interleaved. No need to use it if there are two input files. When two files are given as input and if the input is not correctly interleaved, just print an error message and exit gracefully. If there is a single file input, check while loading, if the sequence names do not look like correctly interleaved, start building the hash table. For FASTQ input the read names are supposed to end with /1 and /2 with the same prefix. For Sanger sequencing data, we may still look at the substrings of the names as FORWARD and REVERSE. 

* Overhaul the global variables, minimize them or remove entirely, by moving them to parameter passes to functions. We may also package them as a struct and pass the pointer to a struct.

* The assembler could be changed, using a config.yaml file. It could be PHRAP for fosmid/plasmid/BAC reads, and SPAdes assembler if the input is short reads from Illumina. The config can also give full path to the assembler binary.

* More automation. The component selection, extraction, and the read fetching from the input file could be built into repeatnet with necessary parameters.

