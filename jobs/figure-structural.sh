#!/usr/bin/env bash

printf "This script compiles the structural plots.\n"
printf "Results from the subplots are manually combined to create the \n"
printf "overall superfigure. "
printf "Consider running the Julia code interactively to inspect the statistics.\n\n"

JOB_PATH=$(pwd)
cd ".."
BASE_PATH=$(pwd)

## Cumulative mass flow plots --------------------------------------------------
DIR="figures/structural/subpanel-a-b/code/"
julia --project=$BASE_PATH -t 1 $BASE_PATH/$DIR/classification.jl
cd $DIR/..
pdflatex -interaction="batchmode" figure-barplot-pathed-versus-looped-efm-flux-hepg2.tex
rm figure-barplot-pathed-versus-looped-efm-flux-hepg2.aux 
rm figure-barplot-pathed-versus-looped-efm-flux-hepg2.log
pdflatex -interaction="batchmode" figure-barplot-tree-classification-carbon.tex
rm figure-barplot-tree-classification-carbon.aux
rm figure-barplot-tree-classification-carbon.log
pdflatex -interaction="batchmode" figure-barplot-tree-classification-carbon-scaled.tex
rm figure-barplot-tree-classification-carbon-scaled.aux
rm figure-barplot-tree-classification-carbon-scaled.log
pdflatex -interaction="batchmode" figure-barplot-tree-classification-nitrogen.tex
rm figure-barplot-tree-classification-nitrogen.aux
rm figure-barplot-tree-classification-nitrogen.log
pdflatex -interaction="batchmode" figure-barplot-tree-classification-nitrogen-scaled.tex
rm figure-barplot-tree-classification-nitrogen-scaled.aux
rm figure-barplot-tree-classification-nitrogen-scaled.log
cd $BASE_PATH

DIR="figures/structural/subpanel-c-d/code/"
julia --project=$BASE_PATH -t 1 $BASE_PATH/$DIR/histogram-curves.jl
cd $DIR/..
pdflatex -interaction="batchmode" figure-efm-length-curves-carbon.tex
rm figure-efm-length-curves-carbon.aux figure-efm-length-curves-carbon.log
pdflatex -interaction="batchmode" figure-efm-length-curves-nitrogen.tex
rm figure-efm-length-curves-nitrogen.aux figure-efm-length-curves-nitrogen.log
cd $BASE_PATH

DIR="figures/structural/subpanel-e-f/code/"
julia --project=$BASE_PATH -t 1 $BASE_PATH/$DIR/jitter.jl
cd $DIR/..
pdflatex -interaction="batchmode" figure-jitter-number-efms-carbon.tex
rm figure-jitter-number-efms-carbon.aux figure-jitter-number-efms-carbon.log
pdflatex -interaction="batchmode" figure-jitter-number-efms-nitrogen.tex
rm figure-jitter-number-efms-nitrogen.aux figure-jitter-number-efms-nitrogen.log
cd $DIR/..

DIR="figures/structural/"
cd $DIR
pdflatex -interaction="batchmode" figure-structural.tex
rm figure-structural.aux figure-structural.log
# ------------------------------------------------------------------------------

printf "SCRIPT COMPLETED.\n\n"

