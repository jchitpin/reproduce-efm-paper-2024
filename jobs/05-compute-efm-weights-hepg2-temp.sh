#!/usr/bin/env bash
THREAD=16

printf "This script computes the carbon and nitrogen EFM weights for the \n"
printf "hepg2 network. Bash variable THREAD is set to $THREAD. Set the \n"
printf "thread count higher if you'd like to reproduce the results faster.\n\n"

JOB_PATH=$(pwd)
cd ".."
BASE_PATH=$(pwd)

julia --project=$BASE_PATH -t $THREAD $BASE_PATH/src/gems/hepg2/03-efm-weights/01-compute-efm-weights-glucose-valine.jl
printf "Enumerated and computed EFM weights for model: hepg2.\n\n"

printf "\n"
printf "WEIGHT IDENTIFICATION SCRIPT COMPLETED.\n\n"

