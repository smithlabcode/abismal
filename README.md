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

(1) make sure you clone with the ``--recursive`` flag, which also
clones the `smithlab_cpp` subdirectory

```
$ cd /where/you_want/the_code
$ git clone --recursive git@github.com:smithlabcode/smithlab_cpp.git
```

(2) Build the `abismal` and `abismalidx` programs:
```
$ make all
$ make install
```

### Installation from a release download

(1) Extract the compressed file to a directory (e.g.
`/path/to/abismal`). Once in the directory, run

```
$ ./configure --prefix=/where/you/want/abismal
$ make all
$ make install
```

This will install `abismal` and `abismalidx` inside the `bin`
directory in the output location.

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

|option|long version |arg type |default|description                                        |
|:-----|:------------|:--------|------:|:--------------------------------------------------|
| -i   | -index      | string  |                   | genome index file [required]          |
| -o   | -outfile    | string  | stdout            | output SAM file                       |
| -m   | -mapstats   | string  | [output].mapstats | mapping statistics output file        |
| -t   | -threads    | integer | 1                 | number of mapping threads             |
| -b   | -batch      | integer | 20,000            | number of reads to load at once       |
| -c   | -candidates | integer | 0                 | maximum candidates for comparison     |
| -p   | -max-mates  | integer | 20                | max number of candidates for mating   |
| -l   | -min-frag   | integer | 32                | run abismal on max sensitivity mode   |
| -L   | -max-frag   | integer | 3,000             | run abismal on max sensitivity mode   |
| -M   | -max-error  | double  | 0.1               | max relative number of errors         |
| -s   | -sensitive  |         |                   | run abismal on max sensitivity mode   |
| -a   | -ambig      |         |                   | report a position for ambiguous reads |
| -P   | -pbat       |         |                   | input follows the PBAT protocol       |
| -R   | -random-pbat|         |                   | input follows the random PBAT protocol|
| -A   | -a-rich     |         |                   | reads are A-rich (SE mode)            |
| -V   | -verbose    |         |                   | print more run info                   |

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

Mapping results are reported in SAM format. Some choices in the output
are explicitly highlighted below:
 * Reads are output identically to how they were read, regardless of
   mapped strand
 * the `NM` tag reports the edit distance between the read and the
   output, specifically the sum of mismatches, insertions and
   deletions to the best mapping position.
 * The `CV` tag reports the assumed bisulfite base used to map the
   read. Reads mapped as A-rich will be reported with `CV:A:A`, and
   reads mapped as T-rich will be reported with `CV:A:T`. This tag is
   independent of the strand the read was mapped to.

### Contacts ###

***Andrew D Smith*** *andrewds@usc.edu*

### Copyright ###

Copyright (C) 2018-2020 Andrew D. Smith and Guilherme de Sena Brandine

Authors: Andrew D. Smith and Guilherme de Sena Brandine

ABISMAL is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your
option) any later version.

ABISMAL is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.
