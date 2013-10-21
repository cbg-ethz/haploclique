#HaploClique
We present a novel method to reconstruct the structure of a viral quasispecies from NGS data.
Our approach can be used to:
 - reconstruct local error-corrected haplotypes and estimate their abundance
 - assemble full-length viral haplotypes
 - detect large deletions and insertions from paired-end data.

## CURRENT STATUS
<b>This is a work in progress project and the goal is to make it a one click solution. </b>
 - [x] Local reconstruction
 - [x] Global assembly
 - [ ] Support data sets with coverage > 1000x
 - [ ] InDel prediction

## INSTALL
###Dependencies
Download [saf](https://github.com/armintoepfer/seqalfixer/releases/) and export its parent directory as $SAF enviroment variable.  
HaploClique depends on [boost](http://www.boost.org/) and [cmake](http://www.cmake.org/). You can install them with a package manager of your choice. For OSX, we recommend [macports](http://www.macports.org/).

###HaploClique
```bash
git clone https://github.com/armintoepfer/haploclique
cd haploclique
sh install-additional-software.sh
mkdir build
cd build
cmake ..
```
If you want to install HaploClique to a non-standard directory, change it with `cmake -DCMAKE_INSTALL_PREFIX=<prefix-path> ..`
```
make
make install
```

## USAGE
HaploClique takes a BAM alignment as input for error-correction:  
 - `haploclique-local reference.fasta alignment.bam` 

The reconstructed local haplotypes are saved as:  
 - __data_cliques_paired_R1.fastq__
 - __data_cliques_paired_R2.fastq__
 - __data_cliques_single.fastq__

To assemble global haplotypes after local reconstruction, execute this command until the number and length of haplotypes converged:  
 - `haploclique-assembly reference.fasta`

Perform at least 50 iterations, like `for i in {1..50}; do haploclique-assembly reference.fasta; done`

#####Contact:
```
Armin Töpfer
armin.toepfer (at) gmail.com
http://www.bsse.ethz.ch/cbg/people/armintoepfer
```
