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
git submodule init && git submodule update
mkdir build
cd build
cmake ..
make
make install
```

### Contributions
 [Armin Töpfer](http://www.armintoepfer.com)  
 [Tobias Marschall](https://bioinf.mpi-inf.mpg.de/homepage/index.php?&account=marschal)
 Bernhard Lang
 Marcel Meyerheim
 
###Contact
```
Armin Töpfer
armin.toepfer (at) gmail.com
http://www.armintoepfer.com
```
