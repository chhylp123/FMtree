FMtree: A fast locating algorithm of FM-indexes for genomic data
============






Introduction
-------  

Here are the implementations of "FMtree: A fast locating algorithm of FM-indexes for genomic data".
This project consists of three locating algorithms: FMtree, Original_s and Original_v. 
We also provide a program named "preprocess" to help users to preprocess the input text (genome
).






How to use
-------



For FMtree, Original_s and Original_v:

* For each method, first make it by giving the command 'make' to generate binary.

* After that, build the index for the input text. For example, consider a text "human.fasta", its index consists of 
"human.fasta.index", "human.fasta.index.bwt", "human.fasta.index.sa" and "human.fasta.index.occ". `Please note that the input text can only include 
the characters which belong to {a, c, g, t, A, C, G, T}`.

* When searching, the pattern must saved in "patterns.txt". `Like the text, patterns cannot includes any characters which 
does not belong to {a, c, g, t, A, C, G, T}`. If some patterns consist of such characters, the results of FMtree will be incorrect. Users can randomly generate patterns using our programs (FMtree/Original_v/Original_s), or generate by themselves.

For the program named "preprocess":


* This program is used to preprocess the text (genome). By utilizing this program, any character which does not belong {a, c, g, t, A, C, G, T} is randomly converted to one of {a, c, g, t, A, C, G, T}. The output text will only consists of the sequence of genome. Other information, like chrome names and chrome length will be removed.


* Usage: ./preprocess --index input_text (Note that the output text will be saved in input_text.not_N).

* Example: ./preprocess --index human.fasta (the output text is saved in human.fasta.not_N).



Note
-------
* We adopt the SACA-K algorithm [1] to build the suffix array, and build BWT from suffix array. As such when building the index, the memory requirement of FMtree, Original_s and Original_v are about 5 times larger than that of the input text.


References
-------


[1] Nong G. Practical linear-time O (1)-workspace suffix sorting for constant alphabets[J]. ACM Transactions on Information Systems (TOIS), 2013, 31(3): 15.
