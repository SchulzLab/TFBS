#! /bin/bash --

#
# simple bash tool to delete data points which deceed or exceed a certain
# threshold. By default all lines that exceeds the given threshold in the
# specified field will be deleted and the result will be printed on stdout.
# The input is expected on stdin.
#
#
# USAGE:
#
# obligatory:
# -t NUMBER       threshold for deletion
# -c NUMBER       column that should be tested for (we start indexing at 1)
# 
#
# optional:
# -o FILE         filepath for output
# -l              delete all lines that DECEED the threshold in the spec. column

OPTIND=1

output_file=""
deceed_flag=0
thresh_flag=0
column_flag=0


while getopts ":t:c:o:l" opt; do

    case "$opt" in

        t)
            thresh_flag=1
            threshold=$OPTARG
            ;;
        c)
            column_flag=1
            column=$OPTARG
            ;;
        o)
            output_file=$OPTARG
            ;;
        l)
            deceed_flag=1
            ;;
        \?)
            echo "ERROR: Invalid option -$OPTARG." >&2
            exit 1
            ;;
        :)
            echo "ERROR: Option -$OPTARG requires an argument." >&2
            exit 1
            ;;

    esac
done

if [ $thresh_flag -eq 0 || $column_flag -eq 0 ] then
    echo "ERROR: no threshold (-t) or column specification (-c) provided." >&2
    exit 1
fi

if [ -f $output_file ] then
    if [ $deceed_flag -eq 0 ] then
        awk '$column<=$threshold' >$output_file
    else
        awk '$column>$threshold' >$output_file
    fi
else
    if [ $deceed_flag -eq 0 ] then
        awk '$column<=$threshold'
    else
        awk '$column>$threshold'
    fi
fi
