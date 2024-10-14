#!/usr/bin/env bash

printf "This script compiles the cumulative mass flow plots.\n"
printf "Manually run the code in the figures/cumulative-mass-flow/code \n"
printf "to tweak each subplot x/y-axes. \n\n"

JOB_PATH=$(pwd)
cd ".."
BASE_PATH=$(pwd)

## Cumulative mass flow plots --------------------------------------------------
DIR="figures/cumulative-mass-flow/code/"
julia --project=$BASE_PATH -t 1 $BASE_PATH/$DIR/figure-cumulative-explained-fluxes.jl
cd $BASE_PATH

DIR="figures/cumulative-mass-flow/subpanels/"
cd $DIR
pdflatex -interaction="batchmode" cumulative-source-met-6-legend.tex
pdflatex -interaction="batchmode" cumulative-source-met-6-log.tex
rm cumulative-source-met-6-legend.aux cumulative-source-met-6-legend.log
rm cumulative-source-met-6-log.aux cumulative-source-met-6-log.log

pdflatex -interaction="batchmode" cumulative-source-met-10-legend.tex
pdflatex -interaction="batchmode" cumulative-source-met-10-log.tex
rm cumulative-source-met-10-legend.aux cumulative-source-met-10-legend.log
rm cumulative-source-met-10-log.aux cumulative-source-met-10-log.log

pdflatex -interaction="batchmode" cumulative-source-met-13-legend.tex
pdflatex -interaction="batchmode" cumulative-source-met-13-log.tex
rm cumulative-source-met-13-legend.aux cumulative-source-met-13-legend.log
rm cumulative-source-met-13-log.aux cumulative-source-met-13-log.log

pdflatex -interaction="batchmode" cumulative-source-met-14-legend.tex
pdflatex -interaction="batchmode" cumulative-source-met-14-log.tex
rm cumulative-source-met-14-legend.aux cumulative-source-met-14-legend.log
rm cumulative-source-met-14-log.aux cumulative-source-met-14-log.log

pdflatex -interaction="batchmode" cumulative-source-met-15-legend.tex
pdflatex -interaction="batchmode" cumulative-source-met-15-log.tex
rm cumulative-source-met-15-legend.aux cumulative-source-met-15-legend.log
rm cumulative-source-met-15-log.aux cumulative-source-met-15-log.log

pdflatex -interaction="batchmode" cumulative-source-met-18-legend.tex
pdflatex -interaction="batchmode" cumulative-source-met-18-log.tex
rm cumulative-source-met-18-legend.aux cumulative-source-met-18-legend.log
rm cumulative-source-met-18-log.aux cumulative-source-met-18-log.log

pdflatex -interaction="batchmode" cumulative-source-met-20-legend.tex
pdflatex -interaction="batchmode" cumulative-source-met-20-log.tex
rm cumulative-source-met-20-legend.aux cumulative-source-met-20-legend.log
rm cumulative-source-met-20-log.aux cumulative-source-met-20-log.log

pdflatex -interaction="batchmode" cumulative-source-met-21-legend.tex
pdflatex -interaction="batchmode" cumulative-source-met-21-log.tex
rm cumulative-source-met-21-legend.aux cumulative-source-met-21-legend.log
rm cumulative-source-met-21-log.aux cumulative-source-met-21-log.log

pdflatex -interaction="batchmode" cumulative-source-met-22-legend.tex
pdflatex -interaction="batchmode" cumulative-source-met-22-log.tex
rm cumulative-source-met-22-legend.aux cumulative-source-met-22-legend.log
rm cumulative-source-met-22-log.aux cumulative-source-met-22-log.log

pdflatex -interaction="batchmode" cumulative-source-met-23-legend.tex
pdflatex -interaction="batchmode" cumulative-source-met-23-log.tex
rm cumulative-source-met-23-legend.aux cumulative-source-met-23-legend.log
rm cumulative-source-met-23-log.aux cumulative-source-met-23-log.log

pdflatex -interaction="batchmode" cumulative-source-met-29-legend.tex
pdflatex -interaction="batchmode" cumulative-source-met-29-log.tex
rm cumulative-source-met-29-legend.aux cumulative-source-met-29-legend.log
rm cumulative-source-met-29-log.aux cumulative-source-met-29-log.log

pdflatex -interaction="batchmode" cumulative-source-met-33-legend.tex
pdflatex -interaction="batchmode" cumulative-source-met-33-log.tex
rm cumulative-source-met-33-legend.aux cumulative-source-met-33-legend.log
rm cumulative-source-met-33-log.aux cumulative-source-met-33-log.log

pdflatex -interaction="batchmode" cumulative-source-met-35-legend.tex
pdflatex -interaction="batchmode" cumulative-source-met-35-log.tex
rm cumulative-source-met-35-legend.aux cumulative-source-met-35-legend.log
rm cumulative-source-met-35-log.aux cumulative-source-met-35-log.log

pdflatex -interaction="batchmode" cumulative-source-met-40-legend.tex
pdflatex -interaction="batchmode" cumulative-source-met-40-log.tex
rm cumulative-source-met-40-legend.aux cumulative-source-met-40-legend.log
rm cumulative-source-met-40-log.aux cumulative-source-met-40-log.log

cd $BASE_PATH

DIR="figures/cumulative-mass-flow/"
cd $DIR
pdflatex -interaction="batchmode" figure-cumulative-combined-legend.tex
pdflatex -interaction="batchmode" figure-cumulative-combined-log.tex
cd $BASE_PATH
# ------------------------------------------------------------------------------

printf "SCRIPT COMPLETED.\n\n"

