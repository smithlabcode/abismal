## ABISMAL ##

**A**nother **Bis**ulfite **M**apping **Al**gorithm (ABISMAL) is
a read mapping program for bisulfite sequencing in DNA methylation
studies.

### Requirements ###

Currently ABISMAL requires a C++ compiler that supports the C++11
standard and OpenMP. The default compiler assumed is g++ (comes with
GCC, available on your Linux or OS X machine). The g++ compiler has
supported the C++11 standard since roughly 2012 (GCC 4.7) so this
should not cause any problems. It also requires an OMP library and
headers to be available, which rarely causes problems. ABISMAL also is
capable of reading input files (FASTQ format) that are gzip
compressed.  This requires that the ZLib library is installed on the
system; this is also rarely a problem.

If you have trouble with the `make` part of the installation procedure
described below, please contact me.

### Installation from a clone of the repo ###

(1) Make sure `smithlab_cpp` source is installed from a clone of that
repo:
```
$ cd /where/you_want/the_code
$ git clone git@github.com:smithlabcode/smithlab_cpp.git
$ export SMITHLAB_CPP=`pwd`/smithlab_cpp
```
If you are using this method, then you do not need to build any of the
`smithlab_cpp` code.

(2) Clone the `abismal` source code repo from Github:
```
$ cd /where/you_want/the_code
$ git clone git@github.com:smithlabcode/abismal.git
```

(3) Build the `abismal` and `abismalidx` programs:
```
$ ./configure --enable-hts --prefix=/where/you/want/the/binaries
$ make all
$ make install
```

### Indexing the genome ###

The index can be constructed as follows, based on a genome existing
entirely in a single FASTA format file:
```
$ abismalidx <genome.fa> <index-file>
```

### Bisulfite mapping ###

single-end reads
```
$ abismal [options] -i <index-file> -o <output-file> <reads.fq>
```
paired-end reads
```
$ abismal [options] -i <index-file> -o <output-file> <read_1.fq> <read_2.fq>
```

### ABISMAL Options ###

|option|long version |arg type |default|description                           |
|:-----|:------------|:--------|------:|:-------------------------------------|
| -i   | -index      | string  |       | index files from abismal [reqd]      |
| -o   | -outfile    | string  |       | output file name [reqd]              |
| -t   | -threads    | integer | 1     | number of threads to use             |
| -m   | -mismatches | integer | 6     | max allowed mismatches               |
| -s   | -shifts     | integer | 3     | number of seed shifts                |
| -b   | -batch      | integer | 1M    | reads to load in RAM at once         |
| -c   | -candidates | integer | 3000  | max candidates for full comparison   |
| -p   | -max-mates  | integer | 20    | max candidates as mates (pe mode)    |
| -l   | -min-frag   | integer | 32    | min fragment size (pe mode)          |
| -L   | -max-frag   | integer | 3000  | max fragment size (pe mode)          |
| -a   | -ambig      |         |       | report a posn for ambiguous mappers  |
| -P   | -pbat       |         |       | input data follow the PBAT protocol  |
| -A   | -a-rich     |         |       | indicates reads are a-rich (se mode) |
| -v   | -verbose    |         |       | print more run info                  |

### Examples ###

(1) **Indexing the genome**

To make an index for hg38:
```
$ abismalidx hg38.fa hg38.abismalidx
```
In the process of building the index, the names of chromosomes will be
truncated at the first whitespace character.

(2) **Bisulfite Mapping**

To map reads to human genome hg38:
```
$ abismal -i hg38.abismalidx -o reads.mr reads.fq
```
### Output Format ###

**The Mapped Read (MR) Format**

This format retains less total information than is typically in a BAM
file, but has all the information required for typical downstream DNA
methylation analysis.
* RNAME (chromosome name)
* SPOS (start position, 0-based)
* EPOS (end position, 0-based)
* QNAME (read name)
* MISMATCH (number of mismatches)
* STRAND (forward or reverse strand)
* QSEQ (the original input read)

If paired-end reads are mapped in proper pair, the QNAME is added
"FRAG:" in the beginning of the read name, the STRAND is the strand of
the first mate mapped and QSEQ is merged according to their mapping
positions. The overlap segment of QSEQ is from the mate R1 or mate R2
and it is the one with less number of 'N' in the read
sequence. MISMATCH is the sum of mismatches in the mate R1 and
mismatches in the mate R2. If paired-end reads are not mapped in
proper pair, they are treated as single-end reads. If the `-a` option
is set, currently there is no indication about which reads map
ambiguously.

### Contacts ###

***Andrew D Smith*** *andrewds@usc.edu*

### Copyright ###

Copyright (C) 2018-2019 Andrew D. Smith

Authors: Andrew D. Smith

ABISMAL is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your
option) any later version.

ABISMAL is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.
