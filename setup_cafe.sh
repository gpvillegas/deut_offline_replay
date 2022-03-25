#!/bin/bash


#echo 'Setting up CaFe ENVIRONMENT variables . . .'


echo '-------------------------------------'
echo 'Checking CaFe directories exist . . .'
echo '-------------------------------------'
echo ''

declare -a arr=("ROOTfiles" "REPORT_OUTPUT", "raw", )


for i in "${arr[@]}"
do
    if [[ -e "$i" ]];
    then
	echo "Directory $i exists." 
	echo ''
    else
	echo "Error: Directory $i does not exists."
	echo "Need to create one or provide a symbolic link "
	echo ''
    fi
    
    # or do whatever with individual element of the array
done

exit 0

