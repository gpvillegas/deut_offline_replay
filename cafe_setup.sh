#!/bin/bash




# Check if necessary symlinks exist in hallc_replay_cafe
sym_link=("hcana" 'raw' 'ROOTfiles' 'REPORT_OUTPUT')

for i in "${sym_link[@]}"
do
    if [[ -L "$i" && -e "$i" ]]; then
	echo "\n=> symbolic link to $i exists and is valid: "
	readlink $i
    else
	echo "\n=> symbolic link to $i does NOT exist or is broken"    
    fi
done

#readlink $sym_link
