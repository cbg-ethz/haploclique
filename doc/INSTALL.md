## INSTALL
###Dependencies
HaploClique depends on [boost](http://www.boost.org/) and
[cmake](http://www.cmake.org/).
You can install them with a package manager of your choice.

####Ubuntu:
```
apt-get install libncurses5-dev cmake libboost-all-dev git build-essential zlib1g-dev
```

####macOS:
Please XCode and its command line tools, and with [brew](http://brew.sh/):
```
brew install cmake boost ninja
```

####Windows:
HaploClique has not been tested on Windows.

###Installation routine:
If you want to install HaploClique to a non-standard directory, change it with `cmake -DCMAKE_INSTALL_PREFIX=<prefix-path> ..`
```bash
git clone https://github.com/cbg-ethz/haploclique && cd haploclique
git submodule update --init --remote
mkdir build && cd build
cmake .. && make
```

Want blazing fast builds? Install [ninja](https://ninja-build.org/)
on your machine and execute
```
cmake -GNinja .. && ninja
```

Run tests?
```
ninja check
```