In the RepeatNet codebase, the **hash** mechanism is used as an efficient lookup and pairing indexing system for sequence reads, allowing the program to quickly pair forward and reverse reads under the same clone name without expensive linear searches.

Here is a detailed breakdown of how the hash is defined, configured, and used:

### 1. The Hash Function and Algorithms
The core hash function is defined in [repeatnet.c:L958-L1073](file:///home/calkan/code/bilkentcompgen/repeatnet/repeatnet.c#L958-L1073):
```c
int hash(char *this_name, int nseq, int type);
```
It supports **15 different string hashing algorithms** (selected via the `type` parameter) to convert clone name strings into an index:
* **Types 0–4**: Custom summation and Rabin-Karp style polynomial hashes (e.g., using `pow(7, i)` or `pow(11, i)`).
* **Type 5**: RS Hash
* **Type 6**: JS Hash
* **Type 7**: ELF Hash
* **Type 8**: BKDR Hash
* **Type 9**: SDBM Hash
* **Type 10**: DJB Hash
* **Type 11** (Default): DEK Hash
* **Type 12**: A basic shift-and-XOR hash
* **Type 13**: FNV Hash
* **Type 14**: AP Hash

### 2. Collision Resolution
The hash function uses **open addressing (linear probing)** to handle collisions:
```c
index = sum % nseq;
while (names[index] != NULL && strcmp(names[index], this_name)){
    collisions++;
    index++;
    if (index == nseq)
      index = 0;
}
return index;
```
If a slot in the global `names` array is already occupied by a different name, it probes forward linearly until it finds a slot that is either empty or already holds the matching clone name.

### 3. Key Use Cases

#### A. Pairing Forward and Reverse Read Sequences
When reading input files (especially interleaved format via [readInterleavedFastq](file:///home/calkan/code/bilkentcompgen/repeatnet/repeatnet.c#L781-L830) and single FASTA files via [readSingleFasta](file:///home/calkan/code/bilkentcompgen/repeatnet/repeatnet.c#L496-L660)), the tool needs to group the forward and reverse reads of each clone together:
1. It parses the read headers and extracts the root clone name (using the `end2clone` helper function to strip suffixes like `.FORWARD` or `.REVERSE`).
2. It calculates the clone's hash index using the chosen hash algorithm.
3. It stores/pairs the sequence ends into three parallel global arrays at that specific `index`:
   - `names[index]`: The unique clone base name.
   - `forwards[index]`: The sequence string of the forward read end.
   - `reverses[index]`: The reverse-complemented sequence string of the reverse read end.

#### B. Benchmarking Hash Functions
Because different datasets and name formats can trigger varying collision rates, the command-line interface provides a benchmarking mode via `-t, --hashtest N`. This tests hash functions `3` through `N-1` on the input sequence names and outputs their elapsed time and collision counts so the user can choose the most optimal hash algorithm.

#### C. Retrieving Clone Names by ID/Vertex
During the graph traversal and analysis steps (such as connected component extraction), vertices in the graph correspond to indices in the hash table. When retrieving sequence names for specific graph vertices using `--clones <FILE>`, the hash values act as identifiers to load and match records from the saved `.names` files.
