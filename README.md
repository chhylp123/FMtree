FMtree: A fast locating algorithm of FM-indexes for genomic data
============






Introduction
-------  

Here are the implementations of "FMtree: A fast locating algorithm of FM-indexes for genomic data".
This project consists of three locating algorithms: FMtree, Original_s and Original_v. 
We also provide a program named "preprocess" to help users to preprocess the input text (genome).






How to use
-------



For FMtree, Original_s and Original_v:

* For each method, first make it by giving the command 'make' to generate binary.

* After that, build the index for the input text. For example, consider a text "human.fasta", its index consists of 
"human.fasta.index", "human.fasta.index.bwt", "human.fasta.index.sa" and "human.fasta.index.occ". `Please note that the input text can only include 
the characters which belong to {a, c, g, t, A, C, G, T}`.

* When searching, the pattern must be saved in "patterns.txt". `Like the text, patterns cannot includes any character which 
does not belong to {a, c, g, t, A, C, G, T}`. If some patterns consist of such characters, the results of FMtree would be incorrect. Users can randomly generate patterns using our programs (FMtree/Original_v/Original_s), or generate by themselves.

* Example: If users want to run FMtree, they should first type "./FMtree". Then FMtree will report the following information:

     ![example1](https://github.com/chhylp123/FMtree/raw/master/example1.png) 

    Users can input 1, 2 and 3 to build index, search patterns and generate random patterns, respectively. 


For the program named "preprocess":


* This program is used to preprocess the text (genome). By utilizing this program, any character which does not belong {a, c, g, t, A, C, G, T} is randomly converted to one of {a, c, g, t, A, C, G, T}. The output text will only consists of the sequence of genome. Other information, like chrome names and chrome length will be removed. Note that the input text (genome) must be formated using .fa or .fasta format.


* Usage: ./preprocess --index input_text (note that the output text will be saved in input_text.not_N).

* Example: ./preprocess --index human.fasta (the output text is saved in human.fasta.not_N).


API
-------
To use the FMtree library, you should first include "FMtree/bwt.h" in your project. After that, compile it with " -O3 -mpopcnt" options and link it with "FMtree/saca-k.cpp" and "FMtree/bwt.cpp".

There are several functions you can call when using FMtree:

- **build**: Create FM-index with sampling distance `compress_sa` from (`*input_refer`)[0..`text_length`-1]. The index will be saved in "`filename`.index", "`filename`.index.bwt", "`filename`.index.sa" and "`filename`.index.occ":


```
unsigned int indenpendent_creadte_index(unsigned int text_length, char** input_refer, unsigned int compress_sa, char* filename);
```


- **load** the index from FM-index with prefix `fileName`:
```
unsigned int load_index(char* filename_prefix);
```

- get the result of **count** query of `pattern`[0..length-1].
```
inline unsigned int count(char* pattern, unsigned int length,
	unsigned int* sp, unsigned int* ep, unsigned int* sp1, unsigned int* ep1)
```


- **save** the index to file called fileName:
```
void save(const char *fileName);
```
- **load** the index from file called fileName:
```
void load(const char *fileName);
```
- **free** memory occupied by index:
```
void free();
```
- get the **index size** in bytes (size in memory):
```
unsigned long long getIndexSize();
```
- get the size in bytes of the text used to build the index:
```
unsigned int getTextSize();
```
- get the result of **count** query:
```
unsigned int count(unsigned char *pattern, unsigned int patternLen);
```







Note
-------
* We adopt the SACA-K algorithm [1] to build the suffix array, and build BWT from suffix array. As such when building the index, the memory requirement of FMtree, Original_s and Original_v is about 5 times larger than that of the input text.

* Please note that FMtree, Original_s and Original_v do not share a same index. For each method, users should build its own index.

References
-------


[1] Nong G. Practical linear-time O (1)-workspace suffix sorting for constant alphabets[J]. ACM Transactions on Information Systems (TOIS), 2013, 31(3): 15.
