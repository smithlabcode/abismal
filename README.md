
[![GitHub Downloads](https://img.shields.io/github/downloads/smithlabcode/abismal/total?style=social)](https://github.com/smithlabcode/abismal/releases/latest)
[![Install with Conda](https://anaconda.org/bioconda/abismal/badges/platforms.svg)](https://anaconda.org/bioconda/abismal)
[![Install with Conda](https://img.shields.io/conda/dn/bioconda/abismal?color=red&label=conda%20downloads&style=flat-square)](https://anaconda.org/bioconda/abismal)

## abismal ##

**A**nother **Bis**ulfite **M**apping **Al**gorithm (abismal) is
a read mapping program for bisulfite sequencing in DNA methylation
studies.

Download the latest stable release
[here](https://github.com/smithlabcode/abismal/releases).

See how to [get started](https://github.com/smithlabcode/abismal/blob/master/docs/MANUAL.md#quick-installation) and the
[full program documentation](https://github.com/smithlabcode/abismal/blob/master/docs/MANUAL.md).

### Requirements ###

Currently abismal requires a C++ compiler that supports the C++11
standard and OpenMP. The default compiler assumed is g++ (comes with
GCC, available on your Linux or macOS machine). The g++ compiler has
supported the C++11 standard since roughly 2012 (GCC 4.7) so this
should not cause any problems. It also requires an OMP library and
headers to be available, which rarely causes problems. Instructions to
get HTSlib, for macOS or Linux systems, can be found below.

If you have trouble with the `make` part of the installation procedure
described below, please contact us via e-mail or through a [GitHub
issue](https://github.com/smithlabcode/abismal/issues).

### Documentation ###

The full documentation for abismal can be found
[here](https://github.com/smithlabcode/abismal/blob/master/docs/MANUAL.md). This
explains the use of each parameter in full detail. Below, after
installation instructions, we describe the most common use cases:
indexing a genome and mapping single-end and paired-end reads.

### Installation on Linux ###

These instructions are for building abismal from source, rather than
obtaining it through a package manager like conda.

These instructions assume you have access to `apt` which is installed
on Ubuntu-based and Debian-based distributions. The only difference
for other linux distributions is how you get the dependencies. Likely
all you need is:
```console
$ sudo apt-get install -y libhts-dev
```
If you don't have adminstrator privileges, there are other options.
If you have the `libhts-dev` installed, to build `abismal` the
following should work:
```console
$ wget https://github.com/smithlabcode/abismal/releases/download/v3.2.2/abismal-3.2.2.tar.gz
$ tar -zxvf abismal-3.2.2.tar.gz
$ cd abismal-3.2.2
$ mkdir build && cd build
$ ../configure --prefix=/where/you/want/abismal
$ make
$ make install
```
Be sure that you have permissions to write files to
`/where/you/want/abismal`.  This will install `abismal`, `abismalidx`
and `simreads` inside the `bin` directory of the specified location.

### Installation on macOS ###

The GitHub repo for abismal includes tests that run on macOS 13
(Ventura), and we use the following steps. Although our tests begin
with a "fresh" macOS installation, they have certain tools already
available. In particular, [Homebrew](https://brew.sh) is already
available and possibly some other tools. Homebrew is necessary as the
first step to get the tools and dependencies:
```console
$ brew update
$ brew install gcc
$ brew install htslib gsl
$ brew list --versions gcc
```
At this point, keep the version of `gcc` in mind, because it will probably
be needed below. If you don't already have `abismal` downloaded, the next
step is to download it. Here we will assume you are using a release rather
than a clone. To build from a clone involves at least one more step.
```console
$ wget https://github.com/smithlabcode/abismal/releases/download/v3.2.2/abismal-3.2.2.tar.gz
$ tar -zxvf abismal-3.2.2.tar.gz
$ cd abismal-3.2.2
```
Finally, these steps build the software:
```console
$ mkdir build && build
$ ../configure \
    --prefix=/path/to/install \
    CXX="g++-13" \
    CPPFLAGS="-I$(brew --prefix)/include" \
    LDFLAGS="-L$(brew --prefix)/lib"
$ make
$ make install
```
Notice the `g++-13` in the `../configure` command. This is the version
number referenced above. If you have a different version number (e.g.,
when gcc-14 is the default), you will need to update that number to
correspond to the major version number. Be sure you have permissions
to write to the directory `/path/to/install`.

### How to get the dependencies through conda ###

If you are on linux and do not have adminstrator privileges to get the
dependencies (e.g., HTSlib), you can get them either by building them
directly from source, or through conda. In particular, for obtaining HTSlib
through conda, do the following:
```
$ conda install -c bioconda htslib
```
as explained here at [htslib](https://anaconda.org/bioconda/htslib).
I used conda obtained through miniconda3, which means the default
location for HTSlib to be installed is `~/miniconda3` and then inside
the `lib` and `include` subdirectores. So once this is done, you can
build `abismal` by replacing the `configure` step in the earlier
explanations by
```console
../configure --prefix=/path/to/install \
    CPPFLAGS="-I${HOME}/miniconda3/include" \
    LDFLAGS="-L${HOME}/miniconda3/lib"
```
Note that you can use this approach with both Linux or macOS, but in
the case of macOS you can replace the `LDFLAGS` and `CPPFLAGS` for
conda, but keep the `CXX` variable. Remember not to use tilde (`~`) in
place of the `${HOME}` variable above. It might work, but shouldn't.

### Installation from a clone of the repo ###

This method is likely only useful if you need the most recent update,
and is not recommended for most users. The only difference from the
above explanations for linux and macos is that you will need to clone
the repo, which means you need `git` installed, and you will also need
to build the sources in place and without much reporting in case of any
problems.
```console
$ cd /where/you_want/the_code
$ git clone --recursive git@github.com:smithlabcode/abismal.git
$ cd abismal
$ make
$ make install
```
If you are building from the source in a cloned repo, you will likely
see other ways to accomplish it by examining the files in the root of
the repo.

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

|option|long version     |arg type |default            | description                           |
|:-----|:----------------|:--------|------------------:|:--------------------------------------|
| -i   | -index          | string  |                   | genome index file                     |
| -g   | -genome         | string  |                   | genome file (FASTA)                   |
| -o   | -outfile        | string  |                   | output file (default SAM format)      |
| -s   | -stats          | string  |                   | mapping statistics output file (YAML) |
| -c   | -max-candidates | integer | 100               | max candidates per seed*              |
| -l   | -min-frag       | integer | 32                | minimum fragment length (PE mode)     |
| -L   | -max-frag       | integer | 3000              | maximum fragment length (PE mode)     |
| -m   | -max-distance   | double  | 0.1               | max relative number of errors         |
| -a   | -ambig          | boolean |                   | report a position for ambiguous reads |
| -P   | -pbat           | boolean |                   | input follows the PBAT protocol       |
| -R   | -random-pbat    | boolean |                   | input follows the random PBAT protocol|
| -A   | -a-rich         | boolean |                   | reads are A-rich (SE mode)            |
| -t   | -threads        | integer | 1                 | number of mapping threads             |
| -v   | -verbose        | boolean |                   | print more run info                   |
| -B   | -bam            | boolean | output SAM format | write output in BAM format            |

\* the max candidates parameter controls the amount of "effort" in
mapping. In the "sensitive" step, which aligns reads with smaller
exact match seeds, abismal skips seeds that retrieves more than `c`
candidates. The higher the value of `c`, the more alignments abismal
performs. Note that abismal still aligns reads to every exact match
hit that spans more than half of the read ("specific step"). The
specific step does not change with the value set by `c`.

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

To map reads in BAM format:
```
$ abismal -B -i hg38.abismalidx -o reads.bam reads.fq
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
NAR Genomics and Bioinformatics, 3(4), lqab115.
```

### Copyright ###

Copyright (C) 2018-2023 Andrew D. Smith and Guilherme de Sena Brandine

Authors: Andrew D. Smith and Guilherme de Sena Brandine

abismal is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your
option) any later version.

abismal is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.
