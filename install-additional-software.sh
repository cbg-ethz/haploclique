#! /usr/bin/env sh
git clone https://github.com/lh3/bwa
cd bwa
make
cd ..

for sw in htslib bcftools samtools; do
	git clone --branch=develop git://github.com/samtools/${sw}.git
	cd ${sw}
	make
	cd ..
done
