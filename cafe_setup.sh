#!/bin/bash


Help()
{
    # Display Help
    echo "-------------------------------------------------------"
    echo "This shell script sets up the necessary symbolic" 
    echo "links (or dir.) for the CaFe experiment based on which "
    echo "machine (ifarm, cdaq, local) the user is at."
    echo ""
    echo "Syntax: ./cafe_setup.sh [-h|f]"
    echo ""
    echo "options:"
    echo "h     Print this help display"
    echo "f     Select filesystem in in which to read/write data " 
    echo "      from the CaFe experiment. This option ONLY applies "
    echo "      if you are running this shell script on ifarm. "
    echo "      The options are: test, volatile, work, group"
    echo "      See https://hallcweb.jlab.org/wiki/index.php/CaFe_Disk_Space "
    echo "      for detailed information on each of these filesystems."
    echo ""
    echo "example: ./cafe_setup.sh -f volatile"
    echo "-------------------------------------------------------"    
}


# initialize machine flags to 0
# (depending on where this script gets called, it will turn ON one of these)
ifarm_flg=0
cdaq_flg=0
local_flg=0


# 'deut' is only specific for my lcoal machine (C. Yero)
# change this to match yout local machine host name
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


# define the optional arguments
while getopts ":hf:" option; do
    case $option in
	h) # display Help
            Help	    
	    exit;;       	
	f) # Enter a filesystem name (only appplies for ifarm)
            fsys=$OPTARG;;	      
	\?) # Invalid option
            echo "Error: Invalid option"
	    Help
            exit;;
    esac
done


#--- define tape allocations ---

# where CaFe raw data output to be replayed will be stored (.dat files)
tape_raw_dir="/mss/hallc/c-cafe-2022/raw"

# tape volume for analysis output (simulation or replay output you want to keep long-term)
tape_analysis_out="/mss/hallc/c-cafe-2022/analysis" 


#=================================
# ifarm
# (off-line experiment analysis
# or testing the replay scripts)
#
# =================================
if [[ ifarm_flg -eq 1 ]]; then

    echo "Setting up CaFe experiment symlinks on ifarm for user: "$USER
    
    if [[ $fsys == "volatile" ]]; then	     
	echo 'Setting up symbolic links to volatile filesystem on ifarm . . .'
	base_dir_voli="/volatile/hallc/c-cafe-2022/"	

	echo "Creating dir $base_dir_voli$USER . . ."
	mkdir $base_dir_voli$USER

	echo "Creating symlink to /mss/hallc/c-cafe-2022/raw"
	ln -sf $tape_raw_dir
	
	echo "Creating dir and symlink to $base_dir_voli$USER/REPORT_OUTPUT . . ."
	mkdir $base_dir_voli$USER"/REPORT_OUTPUT"
	ln -sf $base_dir_voli$USER"/REPORT_OUTPUT"
	
	echo "Creating dir and symlink to $base_dir_voli$USER/ROOTfiles . . ."
	mkdir $base_dir_voli$USER"/ROOTfiles"
	ln -sf $base_dir_voli$USER"/ROOTfiles"
	
	
    elif [[ $fsys == "work" ]]; then	     
	echo 'Setting up symbolic links to work filesystem on ifarm . . .'
	base_dir_work="/work/hallc/c-cafe-2022/"

	echo "Creating dir $base_dir_work$USER . . ."
	mkdir $base_dir_work$USER

	echo "Creating symlink to /mss/hallc/c-cafe-2022/raw"
	ln -sf $tape_raw_dir
	
	echo "Creating dir and symlink to $base_dir_work$USER/REPORT_OUTPUT . . ."
	mkdir $base_dir_work$USER"/REPORT_OUTPUT"
	ln -sf $base_dir_work$USER"/REPORT_OUTPUT"
	
	echo "Creating dir and symlink to $base_dir_work$USER/ROOTfiles . . ."
	mkdir $base_dir_work$USER"/ROOTfiles"
	ln -sf $base_dir_work$USER"/ROOTfiles"
	
	
    elif [[ $fsys == "group" ]]; then	     
	echo 'Setting up symbolic links to group filesystem on ifarm . . .'
	base_dir_group="/group/c-cafe-2022/"

	echo "Creating dir $base_dir_group$USER . . ."
	mkdir $base_dir_group$USER

	echo "Creating symlink to /mss/hallc/c-cafe-2022/raw"
	ln -sf $tape_raw_dir
	
	echo "Creating dir and symlink to $base_dir_group$USER/REPORT_OUTPUT . . ."
	mkdir $base_dir_group$USER"/REPORT_OUTPUT"
	ln -sf $base_dir_group$USER"/REPORT_OUTPUT"
	
	echo "Creating dir and symlink to $base_dir_group$USER/ROOTfiles . . ."
	mkdir $base_dir_group$USER"/ROOTfiles"
	ln -sf $base_dir_group$USER"/ROOTfiles"
	
    elif [[ $fsys == "test" || $fsys == "" ]]; then
	echo 'Setting up test symlinks on ifarm . . .'
	base_dir="/lustre19/expphy/volatile/hallc/c-cafe-2022/test_files"
	raw_dir=$base_dir'/raw'
	ROOTfiles_dir=$base_dir'/ROOTfiles'
	REPORT_OUTPUT_dir=$base_dir'/REPORT_OUTPUT'
	ln -sf $raw_dir raw
	#ls -l raw
	ln -sf $ROOTfiles_dir
	#ls -l ROOTfiles
	ln -sf $REPORT_OUTPUT_dir
	#ls -l REPORT_OUTPUT
	
	
    fi
fi


#===============================
# cdaq cluster
# (online experiment analysis)
#===============================

if [[ cdaq_flg -eq 1 ]]; then

    echo "Setting up CaFe experiment symlinks for online-analysis on Hall C cdaq machine"

    base_dir_cdaq="/net/cdaq/cdaql1data/cdaq/hallc-online-cafe2022/"

    echo "Creating symlink to /mss/hallc/c-cafe-2022/raw"
    ln -sf $tape_raw_dir

    echo "Creating dir and symlink to $base_dir_cdaq/REPORT_OUTPUT . . ."
    mkdir $base_dir_cdaq"/REPORT_OUTPUT"
    ln -sf $base_dir_cdaq"/REPORT_OUTPUT"
    
    echo "Creating dir and symlink to $base_dir_cdaq/ROOTfiles . . ."
    mkdir $base_dir_cdaq"/ROOTfiles"
    ln -sf $base_dir_cdaq"/ROOTfiles"

    echo "Creating dir and symlink to $base_dir_cdaq/PDFs . . ."
    mkdir $base_dir_cdaq"/PDFs"
    ln -sf $base_dir_cdaq"/PDFs"
    
fi


#=============================
# local
# (the user local computer)
#=============================

if [[ local_flg -eq 1 ]]; then
    
    # This function checks if necessary dir. exists, else it creates them 
    dir_arr=("raw" "ROOTfiles" "REPORT_OUTPUT")
    	
    echo 'Checking if necessary directories or symlinks exist in local machine (${USER} ${HOSTNAME} ) . . .'
    
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

