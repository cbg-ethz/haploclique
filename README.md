<h1 align="center">
<img src="https://github.com/armintoepfer/haploclique/blob/master/haploclique_logo.png?raw=true" alt="HaploClique"/></h1>
***
<p align="center">Dr. Armin Töpfer, <a href="http://www.armintoepfer.com">armintoepfer.com</a></p>
***
We present a novel method to reconstruct the structure of a viral quasispecies from NGS data.
Our approach can be used to:
 - reconstruct local error-corrected haplotypes and estimate their abundance
 - assemble full-length viral haplotypes
 - detect large deletions and insertions from paired-end data.


######Reference implementation of: 
Armin Töpfer, Tobias Marschall, Rowena A. Bull, Fabio Luciani, Alexander Schönhuth, and Niko Beerenwinkel.  
<b>[Viral quasispecies assembly via maximal clique finding](http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003515)</b>  
<i>PLOS Computational Biology</i>, 2014.

######Overview
 - [Installation](https://github.com/armintoepfer/haploclique/edit/master/README.md#install)
 - [Error correction](https://github.com/armintoepfer/haploclique/edit/master/README.md#error-correction)
 - [Quasispecies assembly](https://github.com/armintoepfer/haploclique/edit/master/README.md#quasispecies-assembly-of-long-range-haplotypes)
 - [Indel prediction](https://github.com/armintoepfer/haploclique/edit/master/README.md#structural-variation)
 - [Contact](https://github.com/armintoepfer/haploclique/edit/master/README.md#contact)

### INSTALL
######Dependencies
Download [ConsensusFixer](https://github.com/armintoepfer/ConsensusFixer/releases) and [InDelFixer](https://github.com/armintoepfer/InDelFixer/releases) and export its parent directory as $SAF enviroment variable.  
HaploClique depends on [boost](http://www.boost.org/), [gnu parallel](http://www.gnu.org/software/parallel/), and [cmake](http://www.cmake.org/). You can install them with a package manager of your choice.

######Ubuntu:  
```
apt-get install libncurses5-dev cmake libboost-all-dev git build-essential zlib1g-dev parallel
```

######OSX 10.8.x:
Please XCode and its command line tools, and with [macports](http://www.macports.org/):
```
port install cmake boost parallel
```

######Windows:
HaploClique has not been tested on Windows. The scripts depend on the bash shell, awk, and sed.  

######Installation routine:
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

### USAGE
Please use HaploClique in an empty directory.  
Please create an <b>index</b> for your alignment `samtools index alignment.bam`

##### Error correction
For error correction, HaploClique takes a BAM alignment and the reference as input:  
`haploclique-assembly -r ../reference.fasta -i ../alignment.bam` 

######All command-line options:
```bash 
USAGE:     haploclique-assembly options...

OPTIONS:
   -h      Show this message
   -H      Show extended help
   -r      Path to the reference genome (required)
   -i      Path to the alignment in BAM format (required)
```

#####Quasispecies assembly of long-range haplotypes 
For quasispecies assembly, HaploClique has to be executed iteratively. Currently, this procedure cannot be automated to handle any input alignment. Manual assembly can be performed by calling the same command until number of haplotypes has converged:
```bash
haploclique-assembly -r ../reference.fasta -i ../alignment.bam
haploclique-assembly -r ../reference.fasta -i ../alignment.bam
haploclique-assembly -r ../reference.fasta -i ../alignment.bam
```

#####Structural variation
For the prediction of large insertion and deletions, please only use alignments in which *all reads are paired-end* and they are *not allowed to overlap*:  
`haploclique-assembly -r ../reference.fasta -i ../alignment.bam -l`  
The detected indels are saved as indels.vcf.

###Ultra-deep next-generation sequencing data workflow
<p align="center">
  <img src="https://github.com/armintoepfer/haploclique/blob/master/haploclique_workflow.png?raw=true" alt="HaploClique workflow"/>
</p>

### Contributions
 [Armin Töpfer](http://www.armintoepfer.com)  
 [Tobias Marschall](http://homepages.cwi.nl/~tm/)
 
###Contact
```
Armin Töpfer
armin.toepfer (at) gmail.com
http://www.armintoepfer.com
```
