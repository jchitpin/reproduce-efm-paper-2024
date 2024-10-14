#!/usr/bin/env bash

printf "This script compiles the benchmarking plots.\n"
printf "Results from the subplots are manually combined to create the \n"
printf "overall superfigure.\n"
printf "Consider running the Julia code interactively to inspect the statistics.\n\n"

JOB_PATH=$(pwd)
cd ".."
BASE_PATH=$(pwd)

## Cumulative mass flow plots --------------------------------------------------
DIR="figures/benchmarking/code/"
julia --project=$BASE_PATH -t 1 $BASE_PATH/$DIR/scatterplot.jl
julia --project=$BASE_PATH -t 1 $BASE_PATH/$DIR/boxplot-hepg2.jl
cd $BASE_PATH

DIR="figures/benchmarking/subpanel-a/"
cd $DIR
pdflatex -interaction="batchmode" figure-barplot-julia-matlab.tex
rm figure-barplot-julia-matlab.aux figure-barplot-julia-matlab.log
cd $BASE_PATH

DIR="figures/benchmarking/subpanel-b/"
cd $DIR
pdflatex -interaction="batchmode" figure-barplot-julia-matlab-time.tex
rm figure-barplot-julia-matlab-time.aux figure-barplot-julia-matlab-time.log
cd $BASE_PATH

DIR="figures/benchmarking/subpanel-c/"
cd $DIR
pdflatex -interaction="batchmode" figure-scatterplot-julia-benchmark.tex
rm figure-scatterplot-julia-benchmark.aux figure-scatterplot-julia-benchmark.log
pdflatex -interaction="batchmode" figure-scatterplot-julia-benchmark-T2.tex
rm figure-scatterplot-julia-benchmark-T2.aux figure-scatterplot-julia-benchmark-T2.log
pdflatex -interaction="batchmode" figure-scatterplot-julia-benchmark-T4.tex
rm figure-scatterplot-julia-benchmark-T4.aux figure-scatterplot-julia-benchmark-T4.log
pdflatex -interaction="batchmode" figure-scatterplot-julia-benchmark-T8.tex
rm figure-scatterplot-julia-benchmark-T8.aux figure-scatterplot-julia-benchmark-T8.log
pdflatex -interaction="batchmode" figure-scatterplot-julia-benchmark-T16.tex
rm figure-scatterplot-julia-benchmark-T16.aux figure-scatterplot-julia-benchmark-T16.log
pdflatex -interaction="batchmode" figure-scatterplot-julia-benchmark-T32.tex
rm figure-scatterplot-julia-benchmark-T32.aux figure-scatterplot-julia-benchmark-T32.log
pdflatex -interaction="batchmode" figure-scatterplot-julia-benchmark-combined-log-lin.tex
rm figure-scatterplot-julia-benchmark-combined-log-lin.aux figure-scatterplot-julia-benchmark-combined-log-lin.log
pdflatex -interaction="batchmode" figure-scatterplot-julia-benchmark-combined-log-log.tex
rm figure-scatterplot-julia-benchmark-combined-log-log.aux figure-scatterplot-julia-benchmark-combined-log-log.log
cd $BASE_PATH

DIR="figures/benchmarking/subpanel-d/"
cd $DIR
pdflatex -interaction="batchmode" figure-boxplot.tex
rm figure-boxplot.aux figure-boxplot.log
cd $BASE_PATH

DIR="figures/benchmarking/"
cd $DIR
pdflatex -interaction="batchmode" figure-benchmarking.tex
rm figure-benchmarking.aux figure-benchmarking.log
cd $BASE_PATH
# ------------------------------------------------------------------------------

printf "SCRIPT COMPLETED.\n\n"

