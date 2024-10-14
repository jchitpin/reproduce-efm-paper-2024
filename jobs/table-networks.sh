#!/usr/bin/env bash

printf "This script compiles the dataset summary statistics table.\n\n"

JOB_PATH=$(pwd)
cd ".."
BASE_PATH=$(pwd)

## DATASET SUMMARY STATISTICS TABLE --------------------------------------------
julia --project=$BASE_PATH -t 1 $BASE_PATH/tables/table-gem-chmc-summary/table-gem-chmc-summary.jl
julia --project=$BASE_PATH -t 1 $BASE_PATH/tables/table-gem-chmc-summary/miscellaneous-table-statistics.jl # run interactively to view the miscellaneous statistics

cd tables
pdflatex -interaction="batchmode" typeset-tables.tex
rm typeset-tables.aux typeset-tables.log
# ------------------------------------------------------------------------------

printf "SCRIPT COMPLETED.\n\n"

