#!/usr/bin/env bash

printf "This script enumerates the carbon and nitrogen EFMs for the following\n"
printf "genome-scale models (GEMs) using 2, 4, 8, 16, and 32 threads.\n\n"

TMAX=32
THREADS=$(nproc --all)
if [ "$THREADS" -lt "$TMAX" ]; then
    printf "Your CPU does not have 32 threads. Please edit this shell script\n"
    printf "to run this benchmark with a smaller number of threads.\n"
    exit
fi

JOB_PATH=$(pwd)
cd ".."
BASE_PATH=$(pwd)
ADD_PATH="/src/gems/"
CONSTRUCTION="01-dataset-construction"
ENUMERATION="02-efm-enumeration"
FNAME1="01-enumerate-atomic-efms-T2.jl"
FNAME2="01-enumerate-atomic-efms-T4.jl"
FNAME3="01-enumerate-atomic-efms-T8.jl"
FNAME4="01-enumerate-atomic-efms-T16.jl"
FNAME5="01-enumerate-atomic-efms-T32.jl"
# ------------------------------------------------------------------------------

## Hepg2 -----------------------------------------------------------------------
SET="hepg2"
# Enumerate EFMs
julia --project=$BASE_PATH -t 2  $BASE_PATH/$ADD_PATH/$SET/$ENUMERATION/$FNAME1
julia --project=$BASE_PATH -t 4  $BASE_PATH/$ADD_PATH/$SET/$ENUMERATION/$FNAME2
julia --project=$BASE_PATH -t 8  $BASE_PATH/$ADD_PATH/$SET/$ENUMERATION/$FNAME3
julia --project=$BASE_PATH -t 16 $BASE_PATH/$ADD_PATH/$SET/$ENUMERATION/$FNAME4
julia --project=$BASE_PATH -t 32 $BASE_PATH/$ADD_PATH/$SET/$ENUMERATION/$FNAME5
printf "Enumerated EFMs for model: hepg2.\n\n"
# ------------------------------------------------------------------------------

printf "\n"
printf "ENUMERATION BENCHMARKING SCRIPT COMPLETED.\n\n"

