#!/bin/bash


Help()
{
    # Display Help
    echo "-----------------------------------------------"
    echo "This shell script sets up the necessary" 
    echo "symbolic links for the CaFe experiment."
    echo ""
    echo "Syntax: ./cafe_setup.sh [-h|f]"
    echo ""
    echo "options:"
    echo "h     Print this help display"
    echo "f     Select filesystem in in which to read/write data " 
    echo "      from the CaFe experiment. This option only applies "
    echo "	if you are running this shell script on ifarm. "
    echo "      The options are: test, volatile, work, group"
    echo "      See https://hallcweb.jlab.org/wiki/index.php/CaFe_Disk_Space "
    echo "      for detailed information on each of these filesystems."
    echo ""
    echo "example: ./cafe_setup.sh -f volatile"
    echo "-----------------------------------------------"    
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



# ifarm
if [[ ifarm_flg -eq 1 ]]; then
   
    if [[ $fsys == "test" ]]; then
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
	
    elif [[ $fsys == "volatile" ]]; then	     
	echo 'Setting up symbolic links to $fsys filesystem on ifarm . . .'
	base_dir_voli="/volatile/hallc/c-cafe-2022/"
	#check if user dir. exists, else create it
	if [ ! -d $base_dir_voli$USER ]; then
	    echo "Creating dir $base_dir_voli$USER . . ."
	    mkdir $base_dir_voli$USER
	    echo "Creating dir $base_dir_voli$USER/REPORT_OUTPUT . . ."
	    mkdir $base_dir_voli$USER"/REPORT_OUTPUT"
	    echo "Creating dir $base_dir_voli$USER/ROOTfiles . . ."
	    mkdir $base_dir_voli$USER"/ROOTfiles"
	fi
	
	
    elif  [[ $fs == "work" ]]; then	     
	echo 'Setting up symbolic links to $fsys filesystem on ifarm . . .'
	base_dir_work="/work/hallc/c-cafe-2022/"
    elif[[ $fs == "group" ]]; then	     
	echo 'Setting up symbolic links to $fsys filesystem on ifarm . . .'
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
