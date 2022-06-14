## abismal ##


**A**nother **Bis**ulfite **M**apping **Al**gorithm (abismal) is
a read mapping program for bisulfite sequencing in DNA methylation
studies.

### Requirements ###

Currently abismal requires a C++ compiler that supports the C++11
standard and OpenMP. The default compiler assumed is g++ (comes with
GCC, available on your Linux or OS X machine). The g++ compiler has
supported the C++11 standard since roughly 2012 (GCC 4.7) so this
should not cause any problems. It also requires an OMP library and
headers to be available, which rarely causes problems. abismal also is
capable of reading input files (FASTQ format) that are gzip
compressed.  This requires that the ZLib library is installed on the
system; this is also rarely a problem.

If you have trouble with the `make` part of the installation procedure
described below, please contact us via e-mail or through a [GitHub
issue](https://github.com/smithlabcode/abismal/issues).

### Documentation ###

The full documentation for abismal can be found
[here](https://github.com/smithlabcode/abismal/blob/master/docs/MANUAL.md). This explains
the use of each parameter in full detail. Below we describe the most common use cases,
specifically installing the software, indexing a genome and mapping single- and paired-end
reads.

### Installation from a clone of the repo ###

(1) make sure you clone with the ``--recursive`` flag, which also
clones the `smithlab_cpp` subdirectory

```
$ cd /where/you_want/the_code
$ git clone --recursive git@github.com:smithlabcode/abismal.git
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

### abismal options ###

|option|long version   |arg type |default             |description                           |
|:-----|:--------------|:--------|------------------:|:--------------------------------------|
| -i   | -index        | string  |                   | genome index file                     |
| -g   | -genome       | string  |                   | genome file (FASTA)                   |
| -o   | -outfile      | string  | stdout            | output file (SAM)                     |
| -s   | -stats        | string  |                   | mapping statistics output file (YAML) |
| -x   | -sensitive    | boolean |                   | run abismal in max sensitivity mode*  |
| -t   | -threads      | integer | 1                 | number of mapping threads             |
| -l   | -min-frag     | integer | 32                | minimum fragment length (PE mode)     |
| -L   | -max-frag     | integer | 3000              | maximum fragment length (PE mode)     |
| -m   | -max-distance | double  | 0.1               | max relative number of errors         |
| -a   | -ambig        | boolean |                   | report a position for ambiguous reads |
| -P   | -pbat         | boolean |                   | input follows the PBAT protocol       |
| -R   | -random-pbat  | boolean |                   | input follows the random PBAT protocol|
| -A   | -a-rich       | boolean |                   | reads are A-rich (SE mode)            |
| -v   | -verbose      | boolean |                   | print more run info                   |

\* in max sensitivity mode, abismal will not skip frequent k-mers when mapping a read. This
makes abismal 4 to 20 times slower, but may increase the number of mapped reads up to 0.5%.
Run abismal in this mode if you only plan on mapping your dataset once and will do lots of
downstream analyses afterwards, or if you are interested in some highly repetitive regions
of the genome.

### Examples ###

(1) **Indexing the genome**

To make an index for hg38:
```
$ abismalidx hg38.fa hg38.abismalidx
```
In the process of building the index, the names of chromosomes will be
truncated at the first whitespace character.

(2) **Bisulfite Mapping**

To map single-end reads in file `reads.fq` to human genome hg38:
```
$ abismal -i hg38.abismalidx -o reads.sam reads.fq
```

To map paired-end reads in files `reads-1.fq` and `reads-2.fq` to human genome hg38:
```
$ abismal -i hg38.abismalidx -o reads.sam reads-1.fq reads-2.fq
```

To map reads in BAM format (requires [samtools](https://www.htslib.org))
```
$ abismal -i hg38.abismalidx reads.fq | samtools view -b >reads.bam
```

To map reads to human genome without requiring a separate index file
(i.e. run both indexing and mapping simultaneously):
```
$ abismal -g hg38.fa -o reads.sam reads.fq
```

Mapping results are reported in SAM format. Some choices in the output
are explicitly highlighted below:
 * Reads are output identically to how they were read, regardless of
   mapped strand.
 * the `NM` tag reports the edit distance between the read and the
   output, specifically the sum of mismatches, insertions and
   deletions to the best mapping position.
 * The `CV` tag reports the assumed bisulfite base used to map the
   read. Reads mapped as A-rich will be reported with `CV:A:A`, and
   reads mapped as T-rich will be reported with `CV:A:T`. This tag is
   independent of the strand the read was mapped to. If reads are not
   mapped in PBAT or random PBAT mode, the first end will always be
   T-rich and the second end will always be A-rich.

### Contacts ###

***Andrew D Smith*** *andrewds@usc.edu*
***Guilherme Sena*** *desenabr@usc.edu*

### Citation ###
The abismal manuscript is available
[here](https://doi.org/10.1093/nargab/lqab115).
If you used `abismal` to analyze your data, please cite us as follows.
```
de Sena Brandine, G., & Smith, A. D. (2021).
Fast and memory-efficient mapping of short bisulfite sequencing reads using a two-letter alphabet.
NAR genomics and bioinformatics, 3(4), lqab115.
```

### Copyright ###

Copyright (C) 2018-2021 Andrew D. Smith and Guilherme de Sena Brandine

Authors: Andrew D. Smith and Guilherme de Sena Brandine

abismal is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your
option) any later version.

abismal is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.


[![DOI](https://zenodo.org/badge/228294373.svg)](https://zenodo.org/badge/latestdoi/228294373)
[![Install with Conda](https://anaconda.org/bioconda/abismal/badges/installer/conda.svg)](https://anaconda.org/bioconda/abismal)
[![Install with Conda](https://anaconda.org/bioconda/abismal/badges/platforms.svg)](https://anaconda.org/bioconda/abismal)
[![Install with Conda](https://anaconda.org/bioconda/abismal/badges/downloads.svg)](https://anaconda.org/bioconda/abismal)
