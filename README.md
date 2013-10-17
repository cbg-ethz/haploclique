#HaploClique
We present a novel method to reconstruct the structure of a viral quasispecies from NGS data.
Our approach can be used to:
 - reconstruct local error-corrected haplotypes and estimate their abundance
 - assemble full-length viral haplotypes
 - detect large deletions and insertions from paired-end data.
 
The source code and further information, how to apply this software package to your data will be available soon.  
This is a work in progress project.

## INSTALL
###Dependencies
Please install and add to your $PATH: [samtools](https://github.com/samtools/samtools/releases/), [bwa](https://github.com/lh3/bwa), and [seqtk](https://github.com/lh3/seqtk). 
Download [saf](https://github.com/armintoepfer/seqalfixer/releases/) and export its parent directory as $SAF enviroment variable.
HaploClique depends on [boost](http://www.boost.org/) and [cmake](http://www.cmake.org/). You can install them with your package manager of your choice. For OSX, we recommend [macports](http://www.macports.org/).

###HaploClique
```bash
git clone https://github.com/armintoepfer/haploclique
cd haploclique
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
HaploClique takes a BAM alignment as input:  
 - `haploclique alignment.bam` for alignments with coverage <1000x  
 - `haploclique-chunk alignment.bam` for ultra-deep sequencing data sets.

The reconstructed local haplotypes are saved as:  
 - __data_cliques_paired_R1.fastq__
 - __data_cliques_paired_R2.fastq__
 - __data_cliques_single.fastq__

The indel predictions are stored in __indel.vcf__

To assemble global haplotypes, execute this command until the number and length of haplotypes converged:  
 - `haploclique-assembly reference_genome.fasta` if reconstructed haplotypes have coverage <1000x  
 - `haploclique-assembly-chunk reference_genome.fasta` if coverage is >=1000x

#####Contact:
```
Armin Töpfer
armin.toepfer (at) gmail.com
http://www.bsse.ethz.ch/cbg/people/armintoepfer
```
