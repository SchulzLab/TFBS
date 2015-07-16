#! /bin/bash --

if [[ "$#" -eq 0 ]]; then

    sort -s -V -k1,1 -k2,2
else
    sort -s -V -k1,1 -k2,2 $1
fi
