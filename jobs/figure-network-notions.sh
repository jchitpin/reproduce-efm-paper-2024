#!/usr/bin/env bash

printf "This script compiles the introduction figure.\n\n"

JOB_PATH=$(pwd)
cd ".."
BASE_PATH=$(pwd)

## NETWORK HAIRBALLS -----------------------------------------------------------
DIR="figures/network-notions-efms-aefms-v2/"
cd $DIR
pdflatex -interaction="batchmode"  figure-network-notions.tex
rm figure-network-notions.aux figure-network-notions.log
cd $BASE_PATH
# ------------------------------------------------------------------------------

printf "SCRIPT COMPLETED.\n\n"

