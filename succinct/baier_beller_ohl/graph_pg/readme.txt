Graphical pan-genome analysis with compressed suffix trees and the Burrows-Wheeler transform
Uwe Baier
Timo Beller
Enno Ohlebusch



Algorithm:
These are the algorithms proposed in
	Graphical pan-genome analysis with compressed suffix trees and the Burrows-Wheeler transform
	Submitted to BIOINF-2015-1242
Please cite the paper above, if you use one of theses algorithms.



Requirements:
	- a modern c++11 ready compiler such as gcc version 4.7
	- the Succinct Data Structure Library (sdsl-lite)
	  (version d533b2600950b4f878cb063ca0cd1bf340c53df4, maybe newer)



Install:
1. Install the SDSL by the commands:
	git clone git@github.com:simongog/sdsl-lite.git sdsl-lite
	cd sdsl-lite
	./install.sh [PATH]

2. Save the path of the sdsl into the variable SDSLLITE e.g. by:
	export SDSLLITE=[PATH to SDSL]

3. Build the executables by the command:
	make



cst_based:
Has theoretical running time O(n).

Example: ./cst_based.x 10 input.fa output
with 10       = k
     input.fa = the input file in FASTA
     output   = the output filename

Will produce the file output file containing the de Bruijn graph.



bwt_based:
Has theoretical running time O(n log sigma) but in practice fast than cst_based.
The space requirement is roughly 1.5n Byte + size of the compressed de Bruijn graph.

Example: ./bwt_based.x input.fa example kfile.txt
with input.fa = the input file in FASTA
     example  = the output filename
     kfile.txt= a file containing the values of k

Will produce two files for each k in kfile.txt:
     example.k*k*.dot
     example.k*k*.start_nodes.txt
Where the *.dot files contain the de Bruijn graph and the *.start_nodes.txt contains
the list of the start nodes.



Notes:
	The programs should compile with commit d533b2600950b4f878cb063ca0cd1bf340c53df4,
	of the SDSL, but may also work with newer versions of the SDSL.



Limitations:
	The input file must be in FASTA format, especially the input may not contain the 0-byte or 1-byte,
	newlines and characters between '>' and the next newline will be removed.
	k must be smaller than the shortest sequence.
