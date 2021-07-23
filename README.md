# CellReports_TooManyPeaks_analysis

This repository contains a collection of analysis scripts used in our Cell Reports manuscript entitled, “TooManyPeaks identifies drug-resistant-specific regulatory elements from single-cell leukemic epigenomes”. The code includes comparisons, purity benchmark, rare population benchmark, random sample generation, and more.

Run on the following machine specifications:

Ubuntu 20.04
512GiB Memory
Intel(R) Xeon(R) CPU E5-2670 v3 @ 2.30GHz
2 physical processors; 24 cores; 48 threads

* Brief file descriptions

cicero.R: Generic running of Cicero for use with clustering procedures.\n
cistopic.R: Generic running of Cistopic for use with clustering procedures.

cistopic_louvain.R: Generic running of Cistopic with Signac UMAP and Louvain clustering for use with clustering procedures.

cusanovich2018.R: Generic running of Cusanovich2018 for use with clustering procedures.
run_episcanpy.py: Generic running of EpiScanpy for use with clustering procedures.
run_episcanpy_paga.py: Generic running of PAGA built on EpiScanpy for use with clustering procedures.
run_apec.py: Generic running of APEC for use with clustering procedures.
signac.R: Generic running of Signac for use with clustering procedures.
snapatac.R: Generic running of snapatac for use with clustering procedures.
get_mat.R: Auxiliary functions to load single cell matrices in R.

generate_random_buenrostro.hs: Generating random population benchmark data sets from the buenrostro data.
generate_random_satpathy.hs: Generating random population benchmark data sets from the satpathy data.
subsample_cellranger_mat.R: Auxiliary functions to subsample a matrix in R.
rare_population_benchmark.hs: Batch script for running the rare population benchmarks.
clustering-contingency: Calculate the contingency matrix for one of the rare benchmark clustering results.
contingency.hs: Batch script for calculating the contingency matrix for measuring rare benchmark accuracy.

purity.hs: Batch script for calculating the clustering (purity) benchmarks of a data set.
diversity.hs: Calculate different measures for clustering accuracy based on purity.

timing.hs: Batch script for timing benchmarks.
