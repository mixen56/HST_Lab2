#!/bin/bash

TMP=./.gen_tmp
SIZE=${1:-1}        # Mb
PORT=${2:-55555}

rm $TMP 2> /dev/null
touch $TMP
SIZE=$(( SIZE * 1024 )) # convert to Kb

# gen size
matrix_size=$(( RANDOM % 5 + 1 ))
echo $matrix_size >> $TMP

size=0
while [ $size -lt $SIZE ]; do
    # check file size
    size=`ls -l --block-size=K $TMP | awk '{print $5}' | tr -d '[[:alpha:]]'`

    # parallel gen 10 matrices
    for k in {1..20}; do
        # gen matrix
        for (( i = 1; i <= matrix_size; i++ )); do
            for (( i = 1; i <= matrix_size; i++ )); do
                line="$RANDOM.$RANDOM $line"
            done
            line=${line% }      # remove last space
            echo $line >> $TMP
            unset line
        done &
    done
done

