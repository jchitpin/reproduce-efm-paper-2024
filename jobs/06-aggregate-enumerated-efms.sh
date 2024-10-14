#!/usr/bin/env bash
THREAD=1

printf "This script aggregates the saved CHMCs into a single carbon and nitrogen JLD2 file.\n"

JOB_PATH=$(pwd)
cd ".."
BASE_PATH=$(pwd)
ADD_PATH="src/gems-aggregated"
FNAME1="01-aggregate-chmcs.jl"

julia --project=$BASE_PATH -t $THREAD $BASE_PATH/$ADD_PATH/$FNAME1

printf "\n"
printf "CARBON AND NITROGEN CHMC AGGREGATION SCRIPT COMPLETED.\n\n"


