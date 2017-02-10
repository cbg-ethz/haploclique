## This script performs an end to end test on haploclique project. The test input data arabis_large_500_01/reads.bam is a simulated dataset generated from Arabis mosaic virus large satellite RNA reference sequence using Simseq simulator. The expected output is presented in arabis_large_500_01/quasispecies.fasta.fasta and is compared to actual output quasispecies.fasta.fasta.

  $ $__HC_EXE --edge_quasi_cutoff_cliques=0.85 --edge_quasi_cutoff_mixed=0.85 --edge_quasi_cutoff_single=0.8 --min_overlap_cliques=0.6 --min_overlap_single=0.5 --no_singletons --significance=4 --noProb0 $TESTDIR/../data/simulation/arabis_large_500_01/reads.bam > /dev/null;
  $ diff quasispecies.fasta.fasta ${TESTDIR}/../data/simulation/arabis_large_500_01/quasispecies.fasta.fasta
