#!/usr/bin/env bash
THREAD=8

printf "This script converts the HepG2 CHMCs into linked lists for \n"
printf "plotting as Sankey diagrams. Bash variable THREAD is set to $THREAD.\n"
printf "Set the thread count higher if you'd like to reproduce the results faster.\n"
printf "Warning that this script takes a while to run!\n\n"

JOB_PATH=$(pwd)
cd ".."
BASE_PATH=$(pwd)
ADD_PATH="src/gems/hepg2/04-network-visualization"
FNAME1="02-linked-list-reactions.jl"
FNAME2="03-linked-lists.jl"

julia --project=$BASE_PATH -t $THREAD $BASE_PATH/$ADD_PATH/$FNAME1
julia --project=$BASE_PATH -t $THREAD $BASE_PATH/$ADD_PATH/$FNAME2

printf "\n"
printf "COMPLETED LINKED-LISTS TO PLOT SANKEY DIAGRAMS FOR HEPG2 MODEL.\n\n"


