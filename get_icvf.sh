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

e=0;
while IFS= read -r line; do
    if [ "$e" -eq 5 ]; then
        voxel_size=$(echo "$line" | awk '{print $1}')
    fi
    
    if [ $(echo "$line" | awk '{print NF}') -gt 5 ]; then
        r=$(echo "$line" | awk '{print $(NF-1)}')
        s=$(echo "$line" | awk '{print $NF}')
        
        if [ $(echo "$s > 0" | bc -l) -eq 1 ]; then
            R=$(echo "scale=6; $r * sqrt(1.01)" | bc -l)
            area_intra=$(echo "scale=6; $area_intra + 3.14159265359 * $R * $R" | bc -l)
        else
            area_intra=$(echo "scale=6; $area_intra + 3.14159265359 * $r * $r" | bc -l)
        fi
    fi
    
    ((e++))
done < "$file"

total_area=$(echo "$voxel_size * $voxel_size" | bc -l)
icvf=$(echo "scale=6; $area_intra / $total_area" | bc -l)

echo "$icvf"