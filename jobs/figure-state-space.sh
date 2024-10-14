#!/usr/bin/env bash

printf "This script compiles the atomic state space plots.\n\n"

JOB_PATH=$(pwd)
cd ".."
BASE_PATH=$(pwd)

## Cumulative mass flow plots --------------------------------------------------
DIR="figures/atomic-state-space/code"
julia --project=$BASE_PATH -t 1 $BASE_PATH/$DIR/jitter.jl
DIR="figures/atomic-state-space/subpanels"
cd $DIR
pdflatex -interaction="batchmode" figure-jitter-fraction-mc-states-carbon.tex
rm figure-jitter-fraction-mc-states-carbon.aux
rm figure-jitter-fraction-mc-states-carbon.log
pdflatex -interaction="batchmode" figure-jitter-fraction-mc-states-nitrogen.tex
rm figure-jitter-fraction-mc-states-nitrogen.aux
rm figure-jitter-fraction-mc-states-nitrogen.log
cd $BASE_PATH
DIR="figures/atomic-state-space"
cd $DIR
pdflatex -interaction="batchmode" figure-atomic-state-space.tex
rm figure-atomic-state-space.aux figure-atomic-state-space.log
# ------------------------------------------------------------------------------

printf "SCRIPT COMPLETED.\n\n"

