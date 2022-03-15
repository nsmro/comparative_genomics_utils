#!/usr/bin/bash

nb_args=$#
 
if [ $nb_args != 2 ]; then
   echo "        bed2chrom.sh
        converts bed files (at least 6 fields, UCSC standard) to chrom format
        Usage:
            chrom2bed.sh bedfile PREFIX> chromfile"
   exit
else
    chromfile=$1
    awk -v prefix=$2 'BEGIN {FS="\t"} {print prefix"\t"$4"\t"$1"\t"$6"\t"$2"\t"$3}' $chromfile
fi
