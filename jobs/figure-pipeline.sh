#!/usr/bin/env bash

printf "This script compiles the pipeline figure.\n\n"

JOB_PATH=$(pwd)
cd ".."
BASE_PATH=$(pwd)

## NETWORK HAIRBALLS -----------------------------------------------------------
DIR="figures/pipeline/"
cd $DIR
pdflatex -interaction="batchmode"  figure-pipeline.tex
rm figure-pipeline.aux figure-pipeline.log
cd $BASE_PATH
# ------------------------------------------------------------------------------

printf "SCRIPT COMPLETED.\n\n"

