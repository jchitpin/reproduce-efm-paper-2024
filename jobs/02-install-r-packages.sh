#!/usr/bin/env bash

JOB_PATH=$(pwd)
cd ".."
BASE_PATH=$(pwd)

Rscript $BASE_PATH/src/install-r-packages.R

printf "\n"
printf "R PACKAGE INSTALLATION SCRIPT COMPLETED.\n\n"

