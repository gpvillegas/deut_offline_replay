#!/bin/bash

# 12/08/21 - Stephen JD Kay - University of Regina
# Script to fill the run list, executed by "master" script after each run

# 05/25/21 - C. Yero 
# Modified and adpted this script for CaFe 2022 Experiment

#------------
# example of how to
# import pandas as pd

# adding header
# headerList = ['id', 'name', 'profession']

# assuming a csv file has been produced already . . .
#------
# runlist.csv:
# run,run_type,target,target_mass,SHMS_P,SHMS_Angle,HMS_P,HMS_Angle,Ps1,Ps4,Ps5,T1_scal_rate,T4_scal_rate,
# 3288,heep,LH2,1.00794,8.5,8.3, . . .
# 3289,MF,LD2,
# 3290,SRC,Be9, . . .
#------
# read_file = pd.read_csv (r'runlist.csv')
# read_file.to_excel (r'runlist.xlsx', index = None, header=True)

# check if csv file exists, if it does not, create a csv, add a header and add the 1st entry of columns
# if csv file exists, then it already has a header (assuming it was created with this code), and simply add rows

#---------------

# Set up paths depending upon location

# Set path depending upon hostname. Change or add more as needed  
# Note, farm paths are only for testing
if [[ "${HOSTNAME}" = *"cdaq"* ]]; then
    REPLAYPATH="/home/cdaq/c-cafe-2022/cafe_online_replay"
elif [[ "${HOSTNAME}" = *"farm"* ]]; then  
    REPLAYPATH="/w/hallc-scshelf2102/c-cafe-2022/$USER/cafe_online_replay"
elif [[ "${HOSTNAME}" = *"deuteron"* ]]; then
    REPLAYPATH="/Users/nuclear/cafe_online_replay"
elif [[ "${HOSTNAME}" = *"physics"* ]]; then
    REPLAYPATH="/Users/deuteron/CaFe-Online/cafe_online_replay"
fi

# Run number and run type should be read in by the "master" script, automating the target would be more difficult, master script prompts for this
RUN_NUM=$1
EVT_NUM=$2
RUNTYPE=$3  # "sample", "prod", 
KINTYPE=$4  # "bcm_calib", "lumi", "optics", "heep_singles", "heep_coin", "MF", "SRC"  (these are for specific cuts to be applied)

# cafe experiment runlist filename
RUNLIST="${REPLAYPATH}/runlist_cafe_2022.csv"
HEADER="Run,heep_singles,Comment"
# cafe standard kinematics file (to get beam, momenta, target info)
#KINFILE="${REPLAYPATH}/DBASE/COIN/standard.kinematics"

# cafe output
REPORTFILE="${REPLAYPATH}/CAFE_OUTPUT/REPORT/cafe_${RUN_TYPE}_report_${RUN_NUM}_${EVT_NUM}.txt" 

# Get information available in standard.kinematics, execute a python script to do this for us
#KINFILE_INFO=`python3 $REPLAYPATH/UTILS_CAFE/online_scripts/kinfile.py ${KINFILE} ${RUN_NUM}` # The output of this python script is just a comma separated string

# Split the string we get to individual variables, easier for printing and use later
#EBeam=`echo ${KINFILE_INFO}     | cut -d ','  -f1`
#tgt_mass=`echo ${KINFILE_INFO}  | cut -d ','  -f2`
#SHMS_Angle=`echo ${KINFILE_INFO}| cut -d ','  -f3` # Cut the string on , delimitter, select field (f) 1, set variable to output of command
#SHMS_P=`echo ${KINFILE_INFO}    | cut -d ','  -f4`
#SHMS_mass=`echo ${KINFILE_INFO} | cut -d ','  -f5`
#HMS_Angle=`echo ${KINFILE_INFO} | cut -d ','  -f6`
#HMS_P=`echo ${KINFILE_INFO}     | cut -d ','  -f7`
#HMS_mass=`echo ${KINFILE_INFO}  | cut -d ','  -f8`

# Get information available in the report file
if [[ -f ${REPORTFILE} ]]; then
    print_gen_run=`python $REPLAYPATH/UTILS_CAFE/online_scripts/reportfile.py ${REPORTFILE} gen_run_info`
    run_num=`echo ${print_gen_run} | cut -d ',' -f1`
fi

echo "========================================================================="
echo "These values autofill into the run list ..."
echo
echo "Run number: $RUN_NUM"
echo "Run type: $RUNTYPE"
echo "========================================================================="

while true; do
    read -p "Do these values all look correct/reasonable ? (Please answer yes or no) " yn
    case $yn in
        [Yy]* ) break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac
done

# Ask user for a comment
read -p "Enter any comments pertinent to run ${RUN_NUM} (comments will be added to cafe_run_list.csv) : " Comment
Comment=$(echo "$Comment" | tr "," " ") # Remove any commas from the comment line as this will cause... issues
Comment=$(echo "$Comment" | tr ";" " ") # Remove any semicolons from the comment line as well, grammar get out!
Comment=$(echo "$Comment" | tr "\t" " ") # Tabs can go to hell too

# Need to fix widths of entries with blank space at some point, see the test file for widths (based on headers)
RUNLIST_INFO="${RUN_NUM},${RUNTYPE}, $Comment"

# Check if there is already an entry for this run number, if there is, ask if you want to overwrite it, if not, print it to the file
DuplicateLines=() # Array to store line numbers of duplicated entries
LineNum=1 # Counter, starts at 1 since we skip the header

# Read run list, find any lines which already include an entry for this run number
#while IFS='' read -r line || [[ -n "$line" ]]; do
#    LineNum=$(($LineNum + 1 ))
#    if [[ `echo ${line} | cut -d ','  -f1` == ${RUN_NUM} ]]; then
#	DuplicateLines[${#DuplicateLines[@]}]="${LineNum}"
#    fi
#done <  <(tail -n +2 ${RUNLIST}) # Ignores header line by using tail here

echo ${HEADER} >> ${RUNLIST}
echo ${RUNLIST_INFO} >> ${RUNLIST}

#if [[ `echo "${#DuplicateLines[@]}"` != 0 ]]; then # If the array is not empty, check some stuff and check with the user if they want to remove duplicate entries
    # Ask if the user wants to remove duplicate lines, in a grammatically correct manner :)
#    if [[ `echo "${#DuplicateLines[@]}"` == 1 ]]; then 
#	read -p "$(echo "${#DuplicateLines[@]}") entry already found in the runlist for run ${RUN_NUM}, delete dupliacte entry and print new entry to file? <Y/N> " prompt
#    elif [[ `echo "${#DuplicateLines[@]}"` -gt 1 ]]; then
#	read -p "$(echo "${#DuplicateLines[@]}") entries already found in the runlist for run ${RUN_NUM}, delete dupliacte entries and print new entry to file? <Y/N> " prompt
#    fi
#    if [[ $prompt == "y" || $prompt == "Y" || $prompt == "yes" || $prompt == "Yes" ]]; then
#	DeletedLines=0 # Counter to check how many lines we have deleted already
	# Loop over all line numbers identified earlier as being duplicates, delete them with a sed command
#	for value in "${DuplicateLines[@]}"
#	do
#	    LineNum=$(($value-$DeletedLines)) # We need to account for any lines we delete as we go
#	    sed -i "${LineNum}d" ${RUNLIST}
#	    DeletedLines=$(($DeletedLines + 1))
#	done
#	echo ${RUNLIST_INFO} >> ${RUNLIST} # Print the run list info to the file
#    else echo "Will not remove duplicate entries or print new entry to the file, please edit the runlist manually"
#    fi
#else
#    echo ${RUNLIST_INFO} >> ${RUNLIST} # Print the run list info to the file
#fi
