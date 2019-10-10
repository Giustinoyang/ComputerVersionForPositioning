#!/bin/bash

n=`find ./edge_data -type f | wc -l`

for ((i=1;i<n;i++));
do
    ./edge_boxes_demo edge_data/$i.jpg 100
    
    mkdir ./edge_data/$i
    
    for ((j=1;j<=100;j++));
    do
        mv ./edge_data/$i.$j.jpg ./edge_data/$i
        mv ./edge_data/$i.$j.jpg.txt ./edge_data/$i
    done
done

