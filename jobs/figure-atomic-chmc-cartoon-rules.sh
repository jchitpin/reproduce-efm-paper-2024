#!/usr/bin/env bash

printf "This script compiles the atomic CHMC cartoon rule figure.\n"
printf "Subpanels (f) and (g) are manually tweaked in Inkscape. \n"
printf "Product side colour boxes are manually split into two colours. \n\n"

JOB_PATH=$(pwd)
cd ".."
BASE_PATH=$(pwd)

## ATOMIC CHMC CASES FOR ASSIGNING ATOMIC FLUXES -------------------------------
DIR="figures/atomic-chmc-cartoon-rules"
cd $DIR
pdflatex -interaction="batchmode" header.tex
pdflatex -interaction="batchmode" case-1.tex
pdflatex -interaction="batchmode" case-2.tex
pdflatex -interaction="batchmode" case-3.tex
pdflatex -interaction="batchmode" case-4.tex
pdflatex -interaction="batchmode" case-1-bio-example.tex
pdflatex -interaction="batchmode" case-2-bio-example.tex
pdflatex -interaction="batchmode" case-3-bio-example.tex
pdflatex -interaction="batchmode" figure-atomic-chmc-cartoon-rules.tex
rm header.aux header.log
rm case-1.aux case-1.log
rm case-2.aux case-2.log
rm case-3.aux case-3.log
rm case-4.aux case-4.log
rm case-1-bio-example.aux case-1-bio-example.log
rm case-2-bio-example.aux case-2-bio-example.log
rm case-3-bio-example.aux case-3-bio-example.log
rm figure-atomic-chmc-cartoon-rules.aux figure-atomic-chmc-cartoon-rules.log
cd $BASE_PATH
# ------------------------------------------------------------------------------

printf "SCRIPT COMPLETED.\n\n"

