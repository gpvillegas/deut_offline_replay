#!/bin/bash


declare -a my_link=("ROOTfiles" "REPORT_OUTPUT", "raw" )

for i in "{my_link[@]}"
do
    echo "$i"
    if [ -L ${my_link} ] ; then
	if [ -e ${my_link} ] ; then
	    echo "Good link"
	else
	    echo "Broken link"
	fi
    elif [ -e ${my_link} ] ; then
	echo "Not a link"
    else
	echo "Missing"
    fi
done
