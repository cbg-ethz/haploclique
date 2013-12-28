#HaploClique
We present a novel method to reconstruct the structure of a viral quasispecies from NGS data.
Our approach can be used to:
 - reconstruct local error-corrected haplotypes and estimate their abundance
 - assemble full-length viral haplotypes
 - detect large deletions and insertions from paired-end data.

## CURRENT STATUS
<b>This is a work in progress project. The core algorithm is fully functionality, but some helper-scripts are in still in work.</b>
 - [x] Local reconstruction
 - [ ] Global assembly > 1000x
 - [x] Local reconstruction with coverage > 1000x
 - [x] InDel prediction

## INSTALL
###Dependencies
Download [saf](https://github.com/armintoepfer/seqalfixer/releases/) and [InDelFixer](https://github.com/armintoepfer/InDelFixer/releases) and export its parent directory as $SAF enviroment variable.  
HaploClique depends on [boost](http://www.boost.org/), [gnu parallel](http://www.gnu.org/software/parallel/), and [cmake](http://www.cmake.org/). You can install them with a package manager of your choice.

Ubuntu:  
```
apt-get install libncurses5-dev cmake libboost-all-dev git build-essential zlib1g-dev parallel
```

OSX, please XCode and its command line tools, and with [macports](http://www.macports.org/):
```
port install cmake boost parallel
```

###HaploClique
HaploClique has not been tested on Windows. The scripts depend on the bash shell, awk, and sed.  
If you want to install HaploClique to a non-standard directory, change it with `cmake -DCMAKE_INSTALL_PREFIX=<prefix-path> ..`
```bash
git clone https://github.com/armintoepfer/haploclique
cd haploclique
sh install-additional-software.sh
mkdir build
cd build
cmake ..
make
make install
```

## USAGE
HaploClique takes a BAM alignment as input for error-correction:  
 - `haploclique-assembly -r reference.fasta -a alignment.bam` 

Instructions for quasispecies assembly will follow soon and latest upon publication.

The reconstructed haplotypes are saved as:  
 - __data_cliques_paired_R1.fastq__
 - __data_cliques_paired_R2.fastq__
 - __data_cliques_single.fastq__

#####Contact:
```
Armin Töpfer
armin.toepfer (at) gmail.com
http://www.bsse.ethz.ch/cbg/people/armintoepfer
```
