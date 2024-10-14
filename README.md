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
  TeX files (the values will usually be added as comments in the script).
* The beginning of each script will list what is computed/exported.

## Installation

* Download the repository and install all necessary Julia packages by
  running the following commands in your desired installation directory.

## Setting up the software environment

### Julia

Run the following command in the main repository folder to install all
required Julia packages (with the same versions used in this study).

```julia
julia> Pkg.instantiate()
```

If there are any issues, run the following to manually install the packages:

```console
$ bash jobs/01-install-julia-packages.sh
```

### R

The publication-quality Sankey diagrams are made with the `PantaRhei` and
`grid` packages in R.

```console
$ bash jobs/02-install-r-packages.sh
```

### MATLAB

Download [FluxModeCalculator](https://bioinformaticshome.com/db/tool/FluxModeCalculator) to `/some/directory/flux-mode-calculator`.

Add the following line in `startup.mat` to load this automatically load this
package into MATLAB.

```MATLAB
addpath(genpath('/some/directory/flux-mode-calculator'))
```

## Workflows to reproduce the data

### Enumerating the AEFMs (serial)

`$ bash jobs/03-enumerate-efms.sh`

### Enumerating the AEFMs (parallel)

`$ bash jobs/04-enumerate-efms-parallel.sh`

### Computing AEFM weights

`$ bash jobs/05-compute-efm-weights-hepg2.sh`

### Aggregating AEFMs across the five GEMs

`$ bash jobs/06-aggregate-enumerated-efms.sh`

### Reshaping HepG2 AEFMs for interactively visualizing Sankey diagram

`$ bash jobs/07-linked-list-hepg2.sh`

### Interactively visualizing the Sankey diagrams

Run the following script interactively. The resulting Sankey diagrams will
be visualized with Plotly JavaScript which allows you to rearrange the nodes
by clicking and dragging. Change the variables `m` to plot different source
metabolite indices.

`$ julia --project=. figures/hepg2-weights/subpanels-c-f/code/01-figure-sankey-diagram-interactive.jl`

## Workflows to visualize the results

### Figure 01

`$ bash jobs/figure-pipeline.sh`

### Figure 02

`$ bash jobs/figure-benchmarking.sh`

### Figure 03

`$ bash jobs/figure-structural.sh`

### Figure 04

`$ bash jobs/figure-hepg2-weights.sh`

### Supplementary Figure 01

`$ bash jobs/figure-gem-hairballs.sh`

### Supplementary Figure 02

`$ bash jobs/figure-benchmarking-fluxmodecalculator-parallel.sh`

### Supplementary Figure 03

`$ bash jobs/figure-state-space.sh`

### Supplementary Figure 04

`$ bash jobs/figure-cumulative-mass-flow.sh`

### Supplementary Figure 05

`$ bash jobs/figure-atomic-chmc-cartoon-rules.sh`

### Table 01

`$ bash jobs/table-networks.sh`

Consider running `tables/table-gem-chmc-summary/miscellaneous-table-statistics.jl`
interactively to view the statistics printed to console.

## Miscellaneous

Regenerating the source metabolite chemical structures:

`bash jobs/source-metabolite-structures.sh`

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
