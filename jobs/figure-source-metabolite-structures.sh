#!/usr/bin/env bash

printf "This script compiles the source metabolite structures for the\n"
printf "cumulative mass flow plots.\n\n"

JOB_PATH=$(pwd)
cd ".."
BASE_PATH=$(pwd)

## NETWORK HAIRBALLS -----------------------------------------------------------
DIR="figures/source-metabolite-structures/"
cd $DIR
pdflatex -interaction="batchmode" arginine.tex
rm arginine.aux arginine.log
pdflatex -interaction="batchmode" cystine.tex
rm cystine.aux cystine.log
pdflatex -interaction="batchmode" glucose.tex
rm glucose.aux glucose.log
pdflatex -interaction="batchmode" glutamine.tex
rm glutamine.aux glutamine.log
pdflatex -interaction="batchmode" glutamine-carbon-1-2.tex
rm glutamine-carbon-1-2.aux glutamine-carbon-1-2.log
pdflatex -interaction="batchmode" glutamine-carbon-3-4.tex
rm glutamine-carbon-3-4.aux glutamine-carbon-3-4.log
pdflatex -interaction="batchmode" glutamine-carbon-5.tex
rm glutamine-carbon-5.aux glutamine-carbon-5.log
pdflatex -interaction="batchmode" glycine.tex
rm glycine.aux glycine.log
pdflatex -interaction="batchmode" histidine.tex
rm histidine.aux histidine.log
pdflatex -interaction="batchmode" isoleucine.tex
rm isoleucine.aux isoleucine.log
pdflatex -interaction="batchmode" leucine.tex
rm leucine.aux leucine.log
pdflatex -interaction="batchmode" lysine.tex
rm lysine.aux lysine.log
pdflatex -interaction="batchmode" methionine.tex
rm methionine.aux methionine.log
pdflatex -interaction="batchmode" phenylalanine.tex
rm phenylalanine.aux phenylalanine.log
pdflatex -interaction="batchmode" serine.tex
rm serine.aux serine.log
pdflatex -interaction="batchmode" threonine.tex
rm threonine.aux threonine.log
pdflatex -interaction="batchmode" valine.tex
rm valine.aux valine.log
cd $BASE_PATH
# ------------------------------------------------------------------------------

printf "SCRIPT COMPLETED.\n\n"

