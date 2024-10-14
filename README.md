# Atomic elementary flux modes explain the steady state flow of metabolites in flux networks

This repository contains the scripts necessary to regenerate all figures
and results from the manuscript by J. G. Chitpin and T. J. Perkins. Only
the raw data inputs are stored in this repository due to file size
limitations. The beginning of each script will list what is
computed/exported.

## Requirements

The following is required:

* Julia (minimum version 1.10).
* MATLAB (version X) is optional for reproducing the benchmarking studies.
* A shell to run the scripts/workflows.
* TeX distribution (like TeX Live or MiKTeX) to compile figures.

## Notes on reproducibility

* Benchmarks involving MarkovWeightedEFMs.jl and FluxModeCalculator were
  conducted on a Ryzen 5950X with 64 GB of memory (32 GB recommended
  to run all scripts).
* Certain values are not programmatically exported into the figure/table
  TeX files (benchmark results for FluxModeCalculator).
* The beginning of each script will list what is computed/exported.

## Installation

* Download the repository and install all necessary Julia packages by
  running the following commands in your desired installation directory.

## Setting up the Julia/R/MATLAB environment

## Workflows to reproduce the data

### Enumerating the AEFMs (serial)

### Enumerating the AEFMs (parallel)

### Computing AEFM weights

### Aggregating AEFMs across the five GEMs

### Reshaping HepG2 AEFMs for Sankey diagram

### Constructing the Sankey diagrams

## Workflows to visualize the results

### Figure 01

### Figure 02

### Figure 03

### Figure 04

### Supplementary Figure 01

### Supplementary Figure 02

### Supplementary Figure 03

### Supplementary Figure 04

### Supplementary Figure 05

### Table 01



## Reference

Justin G. Chitpin and Theodore J. Perkins,
*Atomic elementary flux modes explain the steady state flow of metabolites in flux networks*.
biorXiv preprint doi: XX.XXXX/XXXX.XX.XX.XXXXXX

## Acknowledgements

We acknowledge the support of the Natural Sciences and Engineering
Research Council of Canada (NSERC), Discovery grant RGPIN-2019-0660 to
T.J.P. J.G.C. was supported by an NSERC CREATE Matrix Metabolomics
Scholarship and an NSERC Alexander Graham Bell Canada Graduate
Scholarship.
