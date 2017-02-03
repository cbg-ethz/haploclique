## This script performs an end to end test on haploclique project. The test input data reads_HIV-1_50_01.bam is a simulated dataset generated from HIV-1 reference sequence using Simseq simulator. The expected output is presented in quasispecies_HIV-1_50_circleci.fasta.fasta and is compared to actual output quasispecies.fasta.fasta.

  $ haploclique_output_path=${HOME}/${CIRCLE_PROJECT_REPONAME}/test/data/simulation
  $ cd $haploclique_output_path 
  $ haploclique_exe_path=../../../build/src
  $ /home/ubuntu/haploclique/build/src/haploclique --edge_quasi_cutoff_cliques=0.85 --edge_quasi_cutoff_mixed=0.85 --edge_quasi_cutoff_single=0.8 --min_overlap_cliques=0.6 --min_overlap_single=0.5 --no_singletons --significance=4 reads_HIV-1_50_01.bam >/dev/null;
  $ diff quasispecies.fasta.fasta quasispecies_HIV-1_50_weakParam_circleci.fasta.fasta > out.txt;
  $ if [ -s out.txt ]; then echo "Different"; else echo "Same"; fi;
  Same
