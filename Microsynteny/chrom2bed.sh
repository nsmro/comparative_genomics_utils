#!/usr/bin/bash

nb_args=$#
 
if [ $nb_args != 1 ]; then
   echo "        chrom2bed.sh
        converts chrom files to bed6 format
        Usage:
            chrom2bed.sh chromfile > bedfile"
   exit
else
    chromfile=$1
    awk 'BEGIN {FS="\t"} {print $3"\t"$5"\t"$6"\t"$2"\t1\t"$4}' $chromfile
fi
