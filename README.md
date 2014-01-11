#HaploClique
We present a novel method to reconstruct the structure of a viral quasispecies from NGS data.
Our approach can be used to:
 - reconstruct local error-corrected haplotypes and estimate their abundance
 - assemble full-length viral haplotypes
 - detect large deletions and insertions from paired-end data.

## CURRENT STATUS
<b>This is a work in progress project. The core algorithm is fully functionality, but some helper-scripts are in still in work.</b>
 - [x] Local reconstruction
 - [x] Global assembly > 1000x
 - [x] Local reconstruction with coverage > 1000x
 - [x] InDel prediction
 - [ ] Check explizit amplicon mode for PacBio data

## INSTALL
###Dependencies
Download [ConsensusFixer](https://github.com/armintoepfer/ConsensusFixer/releases) and [InDelFixer](https://github.com/armintoepfer/InDelFixer/releases) and export its parent directory as $SAF enviroment variable.  
HaploClique depends on [boost](http://www.boost.org/), [gnu parallel](http://www.gnu.org/software/parallel/), and [cmake](http://www.cmake.org/). You can install them with a package manager of your choice.

Ubuntu:  
```
apt-get install libncurses5-dev cmake libboost-all-dev git build-essential zlib1g-dev parallel
```

OSX 10.8.x, please XCode and its command line tools, and with [macports](http://www.macports.org/):
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
Please use HaploClique in an empty directory.
### Error correction
For error correction, HaploClique takes a BAM alignment and the reference as input:  
`haploclique-assembly -r ../reference.fasta -a ../alignment.bam` 

For ultra-deep next-generation sequencing data sets, please set the minimal overlap to 90%:
`haploclique-assembly -r ../reference.fasta -a ../alignment.bam -o 0.9 -j 0.9` 

All command-line options:
```bash 
USAGE:     haploclique-assembly options...

OPTIONS:
   -h      Show this message
   -r      Path to the reference genome (required)
   -a      Path to the alignment in BAM format (required)
   -q      Edge threshold
   -o      Minimal relative overlap
   -c      Maximal coverage
   -s      Minimal clique size
   -p      PacBio alignment with InDelFixer
   -i      InDelFixer alignment
   -t      Number of iterations
   -u      Shuffle in parallelization
   -e      Continue assembly
   -z      Split size
   -g      No parallelization
   -w      Do not re-use singletons
   -n      Do not use alignment prior
   -l      Only predict indels, no haplotypes
```
### Quasispecies assembly of long-range haplotypes 
For quasispecies assembly, please use the helper script `haploclique-assembly-auto`:
`haploclique-assembly-auto -r ../reference.fasta -i ../alignment.bam`

All command-line options:
```bash
USAGE:     haploclique-assembly-auto options...

OPTIONS:
   -h      Show this message
   -r      Path to the reference genome (required)
   -i      Path to the alignment in BAM format (required)
   -l      Whole genome coverage is below 1000x
   -a      PacBio amplicon mode (Currently non-functional)
```

### Structural variation
For the prediction of large insertion and deletions, please only use alignments in which *all reads are paired-end* and they are *not allowed to overlap*:  
`haploclique-assembly -r ../reference.fasta -a ../alignment.bam -l`
The detected indels are saved as indels.vcf.

## Contributions
 [Armin Töpfer](http://www.bsse.ethz.ch/cbg/people/armintoepfer)  
 [Tobias Marschall](http://homepages.cwi.nl/~tm/)
##Contact:
```
Armin Töpfer
armin.toepfer (at) gmail.com
http://www.bsse.ethz.ch/cbg/people/armintoepfer
```
