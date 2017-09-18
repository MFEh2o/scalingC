#!/bin/bash

for file in "$@"
do
	outfile=$(echo $file | sed 's/grb/nc/') 
	
	cdo -f nc copy $file $outfile
done
