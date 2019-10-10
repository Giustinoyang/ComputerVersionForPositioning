#!/bin/bash
let i=1                               # define an incremental variable
mkdir bak                             # make a backup directory
 
for file in *.jpg                     # *.jpg means all jpg files in current directory
do
    cp ${file} bak
    mv ${file} $i.jpg
    echo -e "${file} renamed as ${i}.jpg"
    let i=i+1
done

rm -r bak
