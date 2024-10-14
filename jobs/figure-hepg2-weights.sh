#!/usr/bin/env bash

printf "This script compiles the HepG2 figures.\n\n"

JOB_PATH=$(pwd)
cd ".."
BASE_PATH=$(pwd)

## SUBPANEL A: CUMULATIVE MASS FLOW PLOTS --------------------------------------
DIR="figures/hepg2-weights/"
julia --project=$BASE_PATH -t 1 $BASE_PATH/$DIR/subpanel-a/code/figure-cumulative-explained-fluxes-glutamine.jl
cd $DIR/subpanel-a
pdflatex -interaction="batchmode" cumulative-source-met-14.tex
rm cumulative-source-met-14.aux cumulative-source-met-14.log
# ------------------------------------------------------------------------------

## SUBPANEL B: NETWORK ---------------------------------------------------------
cd $BASE_PATH
julia --project=$BASE_PATH -t 1 $BASE_PATH/$DIR/subpanel-b/code/figure-gephi-graph.jl
cd $DIR/subpanel-b
pdflatex -interaction="batchmode" network-top-5.tex
pdflatex -interaction="batchmode" network-top-10.tex
rm network-top-5.aux network-top-5.log
rm network-top-10.aux network-top-10.log
# ------------------------------------------------------------------------------

## SUBPANEL C-F: SANKEY --------------------------------------------------------
cd $BASE_PATH
#julia --project=$BASE_PATH -t 1 $BASE_PATH/$DIR/subpanel-c-f/code/01-figure-sankey-diagram-interactive.jl # for generating an interactive Sankey diagram with moveable HTML elements; see figures/hepg2-weights/subpanels-c-f/glutamine/interactive-sankey-diagrams/
julia --project=$BASE_PATH -t 1 $BASE_PATH/$DIR/subpanels-c-f/code/02-figure-sankey-diagram-for-R.jl
Rscript $BASE_PATH/$DIR/subpanels-c-f/code/03-figure-sankey-carbon-1-2.R
Rscript $BASE_PATH/$DIR/subpanels-c-f/code/03-figure-sankey-carbon-1-2-unlabelled.R
Rscript $BASE_PATH/$DIR/subpanels-c-f/code/04-figure-sankey-carbon-3-4.R
Rscript $BASE_PATH/$DIR/subpanels-c-f/code/04-figure-sankey-carbon-3-4-unlabelled.R
Rscript $BASE_PATH/$DIR/subpanels-c-f/code/05-figure-sankey-carbon-5.R
Rscript $BASE_PATH/$DIR/subpanels-c-f/code/05-figure-sankey-carbon-5-unlabelled.R

cd $BASE_PATH
cd $DIR/subpanels-c-f/glutamine/carbon-1-2/03-remove-borders
pdflatex -interaction="batchmode" remove-borders.tex
rm remove-borders.aux remove-borders.log

cd $BASE_PATH
cd $DIR/subpanels-c-f/glutamine/carbon-1-2/04-add-labels
pdflatex -interaction="batchmode" add-labels.tex
rm add-labels.aux add-labels.log

cd $BASE_PATH
cd $DIR/subpanels-c-f/glutamine/carbon-3-4/03-remove-borders
pdflatex -interaction="batchmode" remove-borders.tex
rm remove-borders.aux remove-borders.log

cd $BASE_PATH
cd $DIR/subpanels-c-f/glutamine/carbon-3-4/04-add-labels
pdflatex -interaction="batchmode" add-labels.tex 
rm add-labels.aux add-labels.log

cd $BASE_PATH
cd $DIR/subpanels-c-f/glutamine/carbon-5/03-remove-borders
pdflatex -interaction="batchmode" remove-borders.tex
rm remove-borders.aux remove-borders.log

cd $BASE_PATH
cd $DIR/subpanels-c-f/glutamine/carbon-5/04-add-labels
pdflatex -interaction="batchmode" add-labels.tex
rm add-labels.aux add-labels.log

cd $BASE_PATH
cd $DIR/subpanels-c-f/glutamine/network-mapping
pdflatex -interaction="batchmode" network.tex
rm network.aux network.log
# ------------------------------------------------------------------------------

## ASSEMBLING THE SUPERFIGURE --------------------------------------------------
cd $BASE_PATH
cd $DIR
pdflatex -interaction="batchmode" figure-hepg2-weights.tex
rm figure-hepg2-weights.aux figure-hepg2-weights.log
# ------------------------------------------------------------------------------

printf "SCRIPT COMPLETED.\n\n"

