#!/usr/bin/env bash
THREAD=1

printf "This script constructs Catalyst networks for each of the 5 GEMs and \n"
printf "exports the networks as a Graphviz dot and spreadsheet node/edge table.\n"
printf "These files must be manually imported into Gephi to reproduce the \n"
printf "final graphs and manually tweaked to obtain same/similar looking graphs.\n"
printf "After exporting the graphs, the LaTeX code simply stitches the\n"
printf "subpanels into one big figure.\n\n"

JOB_PATH=$(pwd)
cd ".."
BASE_PATH=$(pwd)

## EXPORT NODE/EDGE TABLES -----------------------------------------------------
julia --project=$BASE_PATH -t $THREAD $BASE_PATH/src/gems/e_coli_core/03-network-visualization/01-export-network-for-gephi.jl
julia --project=$BASE_PATH -t $THREAD $BASE_PATH/src/gems/hepg2/04-network-visualization/01-export-network-for-gephi.jl
julia --project=$BASE_PATH -t $THREAD $BASE_PATH/src/gems/iAB_RBC_283/03-network-visualization/01-export-network-for-gephi.jl
julia --project=$BASE_PATH -t $THREAD $BASE_PATH/src/gems/iIT341/03-network-visualization/01-export-network-for-gephi.jl
julia --project=$BASE_PATH -t $THREAD $BASE_PATH/src/gems/iSB619/03-network-visualization/01-export-network-for-gephi.jl
# ------------------------------------------------------------------------------

## NETWORK HAIRBALLS -----------------------------------------------------------
DIR="figures/gem-hairballs/"
cd $DIR
pdflatex -interaction="batchmode"  figure-gem-hairballs.tex
rm figure-gem-hairballs.aux figure-gem-hairballs.log
cd $BASE_PATH
# ------------------------------------------------------------------------------

printf "SCRIPT COMPLETED.\n\n"

