#! /bin/bash --

awk '{ $4=number[split($4,number,".")]; print }' $1
