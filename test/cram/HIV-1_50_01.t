## This script performs an end to end test on haploclique project. The test input data reads_HIV-1_50_01.bam is a simulated dataset generated from HIV-1 reference sequence using Simseq simulator. The expected output is presented in quasispecies_HIV-1_50.fasta.fasta and is compared to actual output quasispecies.fasta.fasta.

$ $__HC_EXE $TESTDIR/../data/simulation/reads_HIV-1_50_01.bam > /dev/null;
$ diff quasispecies.fasta ${TESTDIR}/../data/simulation/quasispecies_HIV-1_50_circleci.fasta
