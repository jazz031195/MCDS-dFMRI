#!/bin/bash

area_intra=0

while getopts ":f:" opt; do
    case $opt in
        f)
            file=$OPTARG
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
done

if [ -z "$file" ]; then
    echo "File path not provided."
    exit 1
fi

e=0

while IFS= read -r line; do
    if [ "$e" -eq 3 ]; then
        icvf=$(echo "$line" | awk '{print $1}')
    fi
    
    ((e++))
done < "$file"

echo "$icvf"
