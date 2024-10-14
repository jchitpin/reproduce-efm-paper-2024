#!/usr/bin/env bash

printf "This script compiles the supplementary benchmarking plot\n"
printf "showing FluxModeCalculator scaling across threads.\n"
printf "The values are manually taken from the MATLAB code under src/gems/.\n\n"

JOB_PATH=$(pwd)
cd ".."
BASE_PATH=$(pwd)

## Cumulative mass flow plots --------------------------------------------------
DIR="figures/benchmarking-fluxmodecalculator-parallel/"
pdflatex -interaction="batchmode" figure-benchmarking-fluxmodecalculator-parallel.tex
rm figure-benchmarking-fluxmodecalculator-parallel.aux
rm figure-benchmarking-fluxmodecalculator-parallel.log
cd $BASE_PATH
# ------------------------------------------------------------------------------

printf "SCRIPT COMPLETED.\n\n"

