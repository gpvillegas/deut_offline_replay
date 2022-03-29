#!/bin/bash


Help()
{
    # Display Help
    echo "-----------------------------------------------"
    echo "This shell script sets up the necessary" 
    echo "symbolic links for the CaFe experiment."
    echo ""
    echo "Syntax: ./cafe_setup.sh [-s]"
    echo ""
    echo "arguments: "
    echo "arg1: test | Use this argument to point to pre-determined test paths/files" 
    echo "if you just want to test the cafe analysis replay works."
    echo "If no argyments are given, it will default to the production paths/files"
    echo ""
    echo "options:"
    echo "example: ./cafe_setup.sh test"
    echo "-----------------------------------------------"    
}

Help

# initialize machine flags to 0
# (depending on where this script gets called, it will turn ON one of these)
ifarm_flg=0
cdaq_flg=0
local_flg=0


# 'deut' is only specific for my lcoal machine (C. Yero)
if echo $HOSTNAME | grep -q "deut"; then
    local_flg=1
elif echo $HOSTNAME | grep -q "ifarm"; then
    ifarm_flg=1
elif echo $HOSTNAME | grep -q "cdaq"; then
    cdaq_flg=1
else
    echo "Did not find machine HOSTNAME. Open and modify this script to add your machine HOSTNAME"
fi

# print to confirm proper machine HOSTNAME was read
echo 'ifarm_flg='$ifarm_flg
echo 'cdaq_flg='$cdaq_flg
echo 'local_flg='$local_flg


#--- read user inputs ---

fs=$1


# ifarm
if [[ ifarm_flg -eq 1 ]]; then
   
    if [[ $fs == "test" ]]; then
	echo 'Setting up test symlinks on ifarm . . .'
	base_dir="/lustre19/expphy/volatile/hallc/c-cafe-2022/test_files"
	raw_dir=$base_dir'/raw'
	ROOTfiles_dir=$base_dir'/ROOTfiles'
	REPORT_OUTPUT_dir=$base_dir'REPORT_OUTPUT'
	ln -sf $raw_dir raw
	#ls -l raw
	ln -sf $ROOTfiles_dir ROOTfiles
	#ls -l ROOTfiles
	ln -sf $REPORT_OUTPUT_dir REPORT_OUTPUT
	#ls -l REPORT_OUTPUT    
	
    elif [[ $fs == "volatile" ]]; then	     
	echo 'Setting up symlinks to $fs on ifarm . . .'
	base_dir_voli="/volatile/hallc/c-cafe-2022/"
    elif  [[ $fs == "work" ]]; then	     
	echo 'Setting up symlinks to $fs on ifarm . . .'
	base_dir_work="/work/hallc/c-cafe-2022/"
    elif[[ $fs == "group" ]]; then	     
	echo 'Setting up symlinks to $fs on ifarm . . .'
	base_dir_group="/group/c-cafe-2022/"
	
	
	#raw_dir="/mss/hallc/c-cafe-2022/raw"
	#ROOTfiles_dir=""
	
    fi
fi


dir_arr=("raw" "ROOTfiles" "REPORT_OUTPUT")


# local
if [[ local_flg -eq 1 ]]; then
    
    echo 'Checking if necessary directories or symlinks exist in my local machine . . .'

    for i in "${dir_arr[@]}"	     
    do     
	if [[ -L "$i" && -d "$i" ]]; then
	    cmd="ls -l $i"
	    echo "$i is a symlink to a directory and it exists:"
	    eval $cmd 
	elif [[ -d "$i" ]]; then
	    echo "/$i directory exists"	
	else
	    echo "$i symlink is broken or /$i dir does not exist. Creating $i directory now . . ."
	    
	    cmd="mkdir $i"
	    echo $cmd
	    eval $cmd
	    echo "done!"
	fi    
    done
fi
