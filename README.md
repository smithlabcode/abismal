[![GitHub Downloads](https://img.shields.io/github/downloads/smithlabcode/abismal/total?style=social)](https://github.com/smithlabcode/abismal/releases/latest)
[![Install with Conda](https://anaconda.org/bioconda/abismal/badges/platforms.svg)](https://anaconda.org/bioconda/abismal)
[![Install with Conda](https://img.shields.io/conda/dn/bioconda/abismal?color=red&label=conda%20downloads&style=flat-square)](https://anaconda.org/bioconda/abismal)

# Abismal

**A**nother **Bis**ulfite **M**apping **Al**gorithm (abismal) is a read
mapping program for bisulfite sequencing in DNA methylation studies.

* Binaries for [macOS](https://github.com/smithlabcode/abismal/releases/download/v3.3.0/abismal-macOS.tar.gz) and [Linux](https://github.com/smithlabcode/abismal/releases/download/v3.3.0/abismal-Linux.tar.gz).
* Latest stable [release](https://github.com/smithlabcode/abismal/releases/latest).
* [Quickstart](#quickstart)
* [Examples](#examples)
* Full documentation [here](https://github.com/smithlabcode/abismal/blob/master/docs/MANUAL.md) (2025-06-13 needs an update)
* Building from [source](#building-from-source)

## Quickstart

* Download
  [macOS](https://github.com/smithlabcode/abismal/releases/download/v3.3.0/abismal-macOS.tar.gz)
  or
  [Linux](https://github.com/smithlabcode/abismal/releases/download/v3.3.0/abismal-Linux.tar.gz)
  binaries, install through [conda](https://anaconda.org/bioconda/abismal) or
  [source](https://github.com/smithlabcode/abismal/releases/download/v3.3.0/abismal-3.3.0.tar.gz).

* Make an index for your reference genome (assuming human hg38):

  ```console
  ./abismal idx hg38.fa hg38.idx
  ```

* Map reads to the reference genome:

  ```console
  ./abismal map -i hg38.idx -o reads.sam reads_1.fq reads_2.fq
  ```

* See
  [documentation](https://github.com/smithlabcode/abismal/blob/master/docs/MANUAL.md)
  for mapping options and output formats.

## Examples

* Make an index from a reference genome in a single FASTA file (here, `hg38.fa`):

  ```console
  ./abismal idx hg8.fa hg38.idx
  ```

* Make an index using 8 threads:

  ```console
  ./abismal idx -t 8 hg38.fa hg38.idx
  ```

* Map single-end reads to hg38 using an index:

  ```console
  ./abismal map -i hg38.idx -o reads.sam reads.fq
  ```

* Map single-end reads with 64 cores (and get very close to 64x speedup):

  ```console
  ./abismal map -i hg38.idx -o reads.sam -t 64 reads.fq
  ```

* Map paired-end reads in `reads_1.fq` and `reads_2.fq`:

  ```console
  ./abismal map -i hg38.idx -o reads.sam reads_1.fq reads_2.fq
  ```

* Get output in BAM format:

  ```console
  ./abismal map -i hg38.idx -B -o reads.bam reads.fq
  ```

* Get mapping statistics in YAML format:

  ```console
  ./abismal map -i hg38.idx -s reads.stats.yaml -o reads.sam reads.fq
  ```

* Skip making the index (make it on the fly):

  ```console
  ./abismal map -g hg38.fa -o reads.sam reads.fq
  ```

* Map reads from PBAT data:

  ```console
  ./abismal map -i hg38.idx -P -o reads.sam reads.fq
  ```

* Map reads from random PBAT data:

  ```console
  ./abismal map -i hg38.idx -R -o reads.sam reads.fq
  ```

Mapping results are reported in SAM format. Some choices in the output are
explicitly highlighted below:

 * Reads are output identically to how they appear in the input FASTQ files,
   regardless of mapped strand.
 * the `NM` tag reports the edit distance between the read and the output,
   specifically the sum of mismatches, insertions and deletions to the best
   mapping position.
 * The `CV` tag reports the assumed bisulfite base used to map the read. Reads
   mapped as A-rich will be reported with `CV:A:A`, and reads mapped as T-rich
   will be reported with `CV:A:T`. This tag is independent of the strand the
   read was mapped to. If reads are not mapped in PBAT or random PBAT mode,
   the first end will always be T-rich and the second end will always be
   A-rich.

## Building from source

If you are here because the binaries don't work for you, please let us know
and we'll try to fix that.

### Linux

These instructions have been tested for Ubuntu 24.04 and Fedora 41. They will
likely work on most APT and RPM-based distributions in 2025.

* Dependencies:

  Ubuntu/Debian

  ```console
  apt-get update && \
  DEBIAN_FRONTEND=noninteractive \
  apt-get install -y --no-install-recommends \
      ca-certificates \
      g++ \
      make \
      zlib1g-dev \
      libhts-dev \
      automake \
      git \
      wget
  ```

  Fedora/Red Hat

  ```console
  dnf update -y && \
  dnf install -y \
      g++ \
      make \
      zlib-devel \
      libhts-devel \
      wget \
      automake \
      awk \
      git
  ```

  The wget is only needed if for the next step, and the automake and git (and
  awk for Fedora/Red Hat) are only needed for the subsequent step using a
  clone. Your machine likely has these already.

  If you don't have admin privileges on your system, you can use Conda to get
  all the dependencies.  Assuming you don't already have conda installed, this
  will get everything you need:

  ```console
  wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh && \
  sh Miniforge3-Linux-x86_64.sh -bsup ${HOME}/miniforge3 && \
  export PATH=${HOME}/miniforge3/bin:$PATH && \
  conda install -y \
      conda-forge::binutils \
      conda-forge::gxx \
      conda-forge::zlib \
      conda-forge::make \
      conda-forge::automake \
      conda-forge::git \
      bioconda::htslib
  ```

* Build from a source release and install in your home directory:

  ```console
  wget https://github.com/smithlabcode/abismal/releases/download/v3.3.0/abismal-3.3.0.tar.gz
  tar -xf abismal-3.3.0.tar.gz
  cd abismal-3.3.0
  ./configure --prefix=${HOME}
  make
  make install
  ```

* Build from a clone and install in your home directory:

  ```console
  git clone --recursive https://github.com/smithlabcode/abismal.git
  cd abismal
  ./autogen.sh
  ./configure --prefix=${HOME}
  make
  make install
  ```

### macOS

* Dependencies:

  Get a compiler through [xcode](https://developer.apple.com/xcode) or get gcc
  from Homebrew by adding `gcc` to the list below.

  ```console
  brew update && \
  brew install \
      zlib \
      htslib \
      automake \
      git \
      wget
  ```

  The wget is only needed if for the next step, and the automake and git are
  only needed for the subsequent step using a clone. I tested these steps on
  GitHub runners for macOS-15 so you might need additional dependencies.

  Conda as explained above will also work, but you need the conda installer
  script for your macOS system.

* Build from a source release and install in your home directory:

  ```console
  wget https://github.com/smithlabcode/abismal/releases/download/v3.3.0/abismal-3.3.0.tar.gz
  tar -xf abismal-3.3.0.tar.gz
  cd abismal-3.3.0
  ./configure CPPFLAGS="-I$(brew --prefix)/include" LDFLAGS="-L$(brew --prefix)/lib" --prefix=${HOME}
  make
  make install
  ```

* Build from a clone and install in your home directory:

  ```console
  git clone --recursive https://github.com/smithlabcode/abismal.git
  cd abismal
  ./autogen.sh
  ./configure CPPFLAGS="-I$(brew --prefix)/include" LDFLAGS="-L$(brew --prefix)/lib" --prefix=${HOME}
  make
  make install
  ```

If you build from source you can enable a mode for mapping very short reads.
Using such reads is discouraged and rarely helpful. Use the `--help` argument
to the configure script to see how to enable this option.

## Contacts

Andrew D Smith *andrewds@usc.edu*

## Citation

The abismal manuscript is available
[here](https://doi.org/10.1093/nargab/lqab115).  If you used `abismal` to
analyze your data, please cite:

```
de Sena Brandine, G., & Smith, A. D. (2021).
Fast and memory-efficient mapping of short bisulfite sequencing reads using a two-letter alphabet.
NAR Genomics and Bioinformatics, 3(4), lqab115.
```
