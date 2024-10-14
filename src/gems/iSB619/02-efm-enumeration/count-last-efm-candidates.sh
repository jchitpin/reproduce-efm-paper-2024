#!/usr/bin/env bash

x=$(grep -Po "(?<=flux modes:\s)\d+" last-efm-candidates.txt)
y=$(grep -Po "\d+(?=\s\\(after adjacency)" last-efm-candidates.txt | awk '{ SUM += $1} END { print SUM }')
sum=$(($x + $y))

echo "Number of EFMs at second last constraint level: " $x
echo "Number of EFMs at last constraint level:        " $y
echo "Total number of EFMs at last constraint level:  "$sum

