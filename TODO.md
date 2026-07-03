# TODO

* Hash table and the hashtest function seem to be for finding the read pairs if the input is not correctly interleaved. I don't remember this well since I coded it 15 years ago, so be careful and check. We can get rid of this and keep the read pairs without the need for a hash table. If the input is not correctly interleaved, just print an error message and exit gracefully. We use the hash values to fetch the connected components later, so we may need to keep the hashing.

* Overhaul the global variables, minimize them or remove entirely, by moving them to parameter passes to functions. We may also package them as a struct and pass the pointer to a struct.

* The assembler could be changed, using a config.yaml file. It could be PHRAP for fosmid/plasmid/BAC reads, and SPAdes assembler if the input is short reads from Illumina. The config can also give full path to the assembler binary.

* More automation. The component selection, extraction, and the read fetching from the input file could be built into repeatnet with necessary parameters.
