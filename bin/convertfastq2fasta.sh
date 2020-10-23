#!/bin/bash

file1=$1
file2=$2
out=$3
cat ${file1} | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' | sed 's/ /_/g' | sed '/^>/s/$/_f1/' | awk '/^>/{sub(/^>/,">Seq"++i"_");}1' > ${out}_1.fasta
cat ${file2} | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' | sed 's/ /_/g' | sed '/^>/s/$/_f2/' | awk '/^>/{sub(/^>/,">Seq"++i"_");}1' > ${out}_2.fasta
