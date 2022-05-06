#!/bin/bash


echo"================================================="
echo "For help on usage, type: ./cafe_setup.sh -help"
echo"================================================="
Help()
{
    # Display Help
    echo "-------------------------------------------------------"
    echo "This shell script automatically sets up the necessary symbolic" 
    echo "links (or dir.) for the CaFe experiment based on which "
    echo "machine (ifarm, cdaq, local) the user is at."
    echo ""
    echo "Syntax: ./cafe_setup.sh [ -h | -f ]"
    echo ""
    echo "options:"
    echo "-help    Print this help display"
    echo ""
    echo "For users on IFARM: "
    echo "-f    ONLY use this option if you are running this shell script on ifarm. "
    echo "      This option selects filesystem in in which to read/write data " 
    echo "      from the CaFe experiment. "
    echo ""
    echo "      The arguments are: test, volatile, work, group"
    echo "      test: this option (default if no argument is provided) will set" 
    echo "      up pre-determined raw/ ROOTfiles/ and REPORT_OUTPUT/ symbolic links "
    echo "      for testing the CaFe replay and analysis scripts using existing data."
    echo ""
    echo "      volatile, work, group: these options will set symbolic links to the corresponding filesystem. "
    echo "      You would want to set these options depending on which stage of the analysis you are in "
    echo "      for example, select volatile if you are in the beginning stages of off-line analysis."
    echo ""
    echo "      See https://hallcweb.jlab.org/wiki/index.php/CaFe_Disk_Space " 
    echo "      for detailed information on each of these filesystems."  
    echo "        "
    echo "example: ./cafe_setup.sh -f volatile"
    echo "-------------------------------------------------------"    
}

set_hcana_link()
{
    if [[ -z $HCANALYZER ]]; then	
	echo "Environment variable: $HCANALYZER does NOT exist. "
	echo "Please make sure to do: source setup.sh(csh) in hcana first. " 
    else
	echo "Creating hcana symbolic link now  . . ."
	ln -sf $HCANALYZER"/hcana"
	ln -sf $HCANALYZER"/libHallC.so"
	ln -sf $HCANALYZER"/libHallC.so.0.90.0"
    fi    
}


# initialize machine flags to 0
# (depending on where this script gets called, it will turn ON one of these)
ifarm_flg=0
cdaq_flg=0



# define the optional arguments
while getopts ":h:f:" option; do
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


if echo $HOSTNAME | grep -q "ifarm"; then
    ifarm_flg=1
elif echo $HOSTNAME | grep -q "cdaq"; then
    cdaq_flg=1
fi


#if [[ ifarm_flg==0 &&  cdaq_flg==0 ]]; then
#    echo "***************************************"
#    echo " Did not recognize remote machine. "
#    echo " Please run: ./cafe_setup.sh -help "
#    echo " for help in running this script."
#    echo "***************************************"
#fi

#--- define tape allocations ---

# where CaFe raw data output to be replayed will be stored (.dat files(
# but these are NOT directly accessible, one would have to look for them in cache.
tape_raw_dir="/mss/hallc/c-cafe-2022/raw"

# tape volume for analysis output (simulation or replay output you want to keep long-term)
tape_analysis_out="/mss/hallc/c-cafe-2022/analysis" 

#--- define cache allocations ---
cache_raw_dir="/cache/hallc/c-cafe-2022/raw/"
cache_analysis_out="/cache/hallc/c-cafe-2022/analysis/"

#=================================
# ifarm
# (off-line experiment analysis
# or testing the replay scripts)
#
# =================================
if [[ ifarm_flg -eq 1 ]]; then

    echo "Checking if necessary directories or symlinks exist in remote machine: " ${USER}"@"${HOSTNAME}". . ."

    # setup the symbolic links to hcana
    set_hcana_link      
   
    
    if [[ $fsys == "volatile" ]]; then	     
	echo 'Setting up symbolic links to volatile filesystem on ifarm . . .'
	base_dir_voli="/volatile/hallc/c-cafe-2022/"	

	echo "Creating dir $base_dir_voli$USER . . ."
	mkdir $base_dir_voli$USER

	echo "Creating symlink to /mss/hallc/c-cafe-2022/raw"
	ln -sf $tape_raw_dir tape

	echo "Creating symlink to /cache/hallc/c-cafe-2022/raw/"
	ln -sf $cache_raw_dir cache
		
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

	echo "Creating symlink to /cache/hallc/c-cafe-2022/raw/"
	ln -sf $cache_raw_dir cache
		
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

	echo "Creating symlink to /cache/hallc/c-cafe-2022/raw/"
	ln -sf $cache_raw_dir cache
		
	echo "Creating dir and symlink to $base_dir_group$USER/REPORT_OUTPUT . . ."
	mkdir $base_dir_group$USER"/REPORT_OUTPUT"
	ln -sf $base_dir_group$USER"/REPORT_OUTPUT"
	
	echo "Creating dir and symlink to $base_dir_group$USER/ROOTfiles . . ."
	mkdir $base_dir_group$USER"/ROOTfiles"
	ln -sf $base_dir_group$USER"/ROOTfiles"
	
    elif [[ $fsys == "test" ]]; then
	echo 'Setting up test symlinks on ifarm for testing cafe replay scripts . . .'
	base_dir="/lustre19/expphy/volatile/hallc/c-cafe-2022/"
	raw_dir=$base_dir'test_raw'  # this is read-only for users (since dir/ was created by cyero to put raw test files)

	base_dir_user="${base_dir}test_output_${USER}/"
	ROOTfiles_dir=${base_dir_user}"ROOTfiles"
	REPORT_OUTPUT_dir=$base_dir_user"REPORT_OUTPUT"

	mkdir $base_dir_user
	mkdir $ROOTfiles_dir
	mkdir $REPORT_OUTPUT_dir

	unlink raw
	ln -sf $raw_dir raw	
	ln -sf $ROOTfiles_dir
	ln -sf $REPORT_OUTPUT_dir

    elif [[ -z $fsys ]]; then
	echo "No optional argumnet provided. Will default to setting up the symbolic links "
	echo "  to the user directory in volatile for testing cafe replay scripts . . ."
	echo ""
	echo "----------------------------------------------------------------------"
	echo " For help using additional options, please run: ./cafe_setup.sh -help "
	echo "----------------------------------------------------------------------"

	base_dir="/lustre19/expphy/volatile/hallc/c-cafe-2022/"
	raw_dir=$base_dir'test_raw'  # this is read-only for users (since dir/ was created by cyero to put raw test files)

	base_dir_user="${base_dir}test_output_${USER}/"
	ROOTfiles_dir=${base_dir_user}"ROOTfiles"
	REPORT_OUTPUT_dir=$base_dir_user"REPORT_OUTPUT"

	mkdir $base_dir_user
	mkdir $ROOTfiles_dir
	mkdir $REPORT_OUTPUT_dir

	unlink raw
	ln -sf $raw_dir raw	
	ln -sf $ROOTfiles_dir
	ln -sf $REPORT_OUTPUT_dir

	
    fi
fi


#===============================
# cdaq cluster
# (online experiment analysis)
#===============================

if [[ cdaq_flg -eq 1 ]]; then

    echo "Checking if necessary directories or symlinks exist in remote machine: " ${USER}"@"${HOSTNAME}". . ."

    # setup the symbolic links to hcana
    set_hcana_link
    
    base_dir_cdaq="/net/cdaq/cdaql1data/cdaq/hallc-online-cafe2022"

    echo "Creating symlink to /mss/hallc/c-cafe-2022/raw"
    ln -sf $tape_raw_dir

    echo "Creating symlink to /cache/hallc/c-cafe-2022/raw/"
    ln -sf $cache_raw_dir cache
	
    echo "Creating dir and symlink to $base_dir_cdaq/REPORT_OUTPUT . . ."
    mkdir $base_dir_cdaq"/REPORT_OUTPUT"
    ln -sf $base_dir_cdaq"/REPORT_OUTPUT"
    
    echo "Creating dir and symlink to $base_dir_cdaq/ROOTfiles . . ."
    mkdir $base_dir_cdaq"/ROOTfiles"
    ln -sf $base_dir_cdaq"/ROOTfiles"

    echo "Creating dir and symlink to $base_dir_cdaq/HISTOGRAMS . . ."
    mkdir $base_dir_cdaq"/HISTOGRAMS"
    ln -sf $base_dir_cdaq"/HISTOGRAMS"
    
fi


#=============================
# local
# (the user local computer)
#=============================

# assume user is local if NOT on cdaq or ifarm
if [[ ifarm_flg==0 && cdaq_flg==0 ]]; then

    
    # This function checks if necessary dir. exists, else it creates them 
    dir_arr=("raw" "ROOTfiles" "REPORT_OUTPUT" "HISTOGRAMS")
    	
    echo "Checking if necessary directories or symlinks exist in local machine: " ${USER}"@"${HOSTNAME}". . ."

    # setup the symbolic links to hcana
    set_hcana_link

    
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
