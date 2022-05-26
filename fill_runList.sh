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

#---------------

# Set up paths depending upon location

# Set path depending upon hostname. Change or add more as needed  
# Note, farm paths are only for testing
if [[ "${HOSTNAME}" = *"cdaq"* ]]; then
    REPLAYPATH="/home/cdaq/c-cafe-2022/cafe_online_replay"
if [[ "${HOSTNAME}" = *"farm"* ]]; then  
    REPLAYPATH="/w/hallc-scshelf2102/c-cafe-2022/$USER/cafe_online_replay"
elif [[ "${HOSTNAME}" = *"deuteron"* ]]; then
    REPLAYPATH="/Users/nuclear/cafe_online_replay"
fi

# Run number and run type should be read in by the "master" script, automating the target would be more difficult, master script prompts for this
RUNNUMBER=$1
RUNTYPE=$2  #"sample", "prod", "lumi", "tgt_boil", "p_abs"
#KINTYPE=$3  $ "heep", "MF", "SRC"  (these are for specific cuts to be applied)
#TARGET=$3
RUNLIST="${REPLAYPATH}/UTILS_CAFE/runlist_cafe_2022.csv"

# cafe standard kinematics file (to get beam, momenta, target info)
KINFILE="${REPLAYPATH}/DBASE/COIN/standard.kinematics"

# Get cafe report file based upon run type
REPORTFILE="${REPLAYPATH}/CAFE_OUTPUT/REPORT/cafe_${RUN_TYPE}_report_${RUNNUMBER}_-1.txt" 

# Get information available in standard.kinematics, execute a python script to do this for us
KINFILE_INFO=`python3 $REPLAYPATH/UTILS_CAFE/online_scripts/kinfile.py ${KINFILE} ${RUNNUMBER}` # The output of this python script is just a comma separated string

# Split the string we get to individual variables, easier for printing and use later
EBeam=`echo ${KINFILE_INFO}     | cut -d ','  -f1`
tgt_mass=`echo ${KINFILE_INFO}  | cut -d ','  -f2`
SHMS_Angle=`echo ${KINFILE_INFO}| cut -d ','  -f3` # Cut the string on , delimitter, select field (f) 1, set variable to output of command
SHMS_P=`echo ${KINFILE_INFO}    | cut -d ','  -f4`
SHMS_mass=`echo ${KINFILE_INFO} | cut -d ','  -f5`
HMS_Angle=`echo ${KINFILE_INFO} | cut -d ','  -f6`
HMS_P=`echo ${KINFILE_INFO}     | cut -d ','  -f7`
HMS_mass=`echo ${KINFILE_INFO}  | cut -d ','  -f8`

# Get information available in the report file
if [[ -f ${REPORTFILE} ]]; then
    if [[ ${RUNTYPE} != "HeePSing" ]]; then
	REPORTFILE_INFO=`python3 $REPLAYPATH/UTIL_PION/scripts/runlist/reportfile.py ${REPORTFILE}`
    elif [[ ${RUNTYPE} == "HeePSing" ]]; then
	REPORTFILE_INFO=`python3 $REPLAYPATH/UTIL_PION/scripts/runlist/reportfile_HeePSing.py ${REPORTFILE}`
    fi
    Current=`echo ${REPORTFILE_INFO} | cut -d ',' -f1`
    PS1=`echo ${REPORTFILE_INFO} | cut -d ',' -f2`
    PS4=`echo ${REPORTFILE_INFO} | cut -d ',' -f3`
    PS5=`echo ${REPORTFILE_INFO} | cut -d ',' -f4`
    HMS_Rate=`echo ${REPORTFILE_INFO} | cut -d ',' -f5`
    SHMS_Rate=`echo ${REPORTFILE_INFO} | cut -d ',' -f6`
    COIN_Rate=`echo ${REPORTFILE_INFO} | cut -d ',' -f7`
    Charge=`echo ${REPORTFILE_INFO} | cut -d ',' -f8`
    Raw_HMS=`echo ${REPORTFILE_INFO} | cut -d ',' -f9`
    Raw_SHMS=`echo ${REPORTFILE_INFO} | cut -d ',' -f10`
    Raw_COIN=`echo ${REPORTFILE_INFO} | cut -d ',' -f11`
    EDTM=`echo ${REPORTFILE_INFO} | cut -d ',' -f12`
    Tracking=`echo ${REPORTFILE_INFO} | cut -d ',' -f13`
elif [[ ! -f ${REPORTFILE} ]]; then
    Current="ERROR"
    PS1="ERROR"
    PS4="ERROR"
    PS5="ERROR"
    HMS_Rate="ERROR"
    SHMS_Rate="ERROR"
    COIN_Rate="ERROR"
    Charge="ERROR"
    Raw_HMS="ERROR"
    Raw_SHMS="ERROR"
    Raw_COIN="ERROR"
    EDTM="ERROR"
    Tracking="ERROR"
fi

echo "========================================================================="
echo "These values autofill into the run list ..."
echo
echo "Run number: $RUNNUMBER"
echo "Run type: $RUNTYPE"
echo "Target: $TARGET"
echo "Beam energy: $EBeam"
echo "SHMS momentum: $SHMS_P"
echo "SHMS angle : $SHMS_Angle"
echo "SHMS particle mass : $SHMS_mass"
echo "HMS momentum: $HMS_P"
echo "HMS angle: $HMS_Angle"
echo "HMS particle mass : $HMS_mass"
echo "Current: $Current"
if [[ ${RUNTYPE} != "HeePSing" ]]; then
    echo "PS1 : $PS1"
else echo "PS2 : $PS1" # For HeepSingles we care about both ELReal triggers (2 and 4)
fi
echo "PS4 : $PS4"
echo "PS5 : $PS5"
echo "HMS rate [kHz]: $HMS_Rate"
echo "SHMS rate [kHz]: $SHMS_Rate"
if [[ ${RUNTYPE} != "HeePSing" ]]; then
    echo "COIN rate [kHz]: $COIN_Rate"
else echo "COIN rate [kHz] (No coin - SHMS rate duplicated): $COIN_Rate"
fi
echo "Charge [mC]: $Charge"
echo "Raw HMS: $Raw_HMS"
echo "Raw SHMS: $Raw_SHMS"
if [[ ${RUNTYPE} != "HeePSing" ]]; then
    echo "Raw coin: $Raw_COIN"
else echo "Raw coin (No coin - SHMS number duplicated): $Raw_COIN"
fi
echo "EDTM: $EDTM"
if [[ ${RUNTYPE} != "HeePSing" ]]; then
    echo "SHMS hadron tracking: $Tracking"
else echo "SHMS electron tracking: $Tracking"
fi
echo "========================================================================="

while true; do
    read -p "Do these values all look correct? (Please answer yes or no) " yn
    case $yn in
        [Yy]* ) break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac
done

# Ask user for a comment
read -p "Enter number of pi/n events and/or any other comments: " Comment
Comment=$(echo "$Comment" | tr "," " ") # Remove any commas from the comment line as this will cause... issues
Comment=$(echo "$Comment" | tr ";" " ") # Remove any semicolons from the comment line as well, grammar get out!
Comment=$(echo "$Comment" | tr "\t" " ") # Tabs can go to hell too
# Need to fix widths of entries with blank space at some point, see the test file for widths (based on headers)
RUNLIST_INFO="${RUNNUMBER},${RUNTYPE},${TARGET},${EBeam},${SHMS_P},${SHMS_Angle},${HMS_P},${HMS_Angle},${Current},${PS1},${PS4},${PS5},${HMS_Rate},${SHMS_Rate},${COIN_Rate},${Charge},${Raw_HMS},${Raw_SHMS},${Raw_COIN},${EDTM},${Tracking},${Comment}"

# Check if there is already an entry for this run number, if there is, ask if you want to overwrite it, if not, print it to the file
DuplicateLines=() # Array to store line numbers of duplicated entries
LineNum=1 # Counter, starts at 1 since we skip the header
# Read run list, find any lines which already include an entry for this run number
while IFS='' read -r line || [[ -n "$line" ]]; do
    LineNum=$(($LineNum + 1 ))
    if [[ `echo ${line} | cut -d ','  -f1` == ${RUNNUMBER} ]]; then
	DuplicateLines[${#DuplicateLines[@]}]="${LineNum}"
    fi
done <  <(tail -n +2 ${RUNLIST}) # Ignores header line by using tail here

if [[ `echo "${#DuplicateLines[@]}"` != 0 ]]; then # If the array is not empty, check some stuff and check with the user if they want to remove duplicate entries
    # Ask if the user wants to remove duplicate lines, in a grammatically correct manner :)
    if [[ `echo "${#DuplicateLines[@]}"` == 1 ]]; then 
	read -p "$(echo "${#DuplicateLines[@]}") entry already found in the runlist for run ${RUNNUMBER}, delete dupliacte entry and print new entry to file? <Y/N> " prompt
    elif [[ `echo "${#DuplicateLines[@]}"` -gt 1 ]]; then
	read -p "$(echo "${#DuplicateLines[@]}") entries already found in the runlist for run ${RUNNUMBER}, delete dupliacte entries and print new entry to file? <Y/N> " prompt
    fi
    if [[ $prompt == "y" || $prompt == "Y" || $prompt == "yes" || $prompt == "Yes" ]]; then
	DeletedLines=0 # Counter to check how many lines we have deleted already
	# Loop over all line numbers identified earlier as being duplicates, delete them with a sed command
	for value in "${DuplicateLines[@]}"
	do
	    LineNum=$(($value-$DeletedLines)) # We need to account for any lines we delete as we go
	    sed -i "${LineNum}d" ${RUNLIST}
	    DeletedLines=$(($DeletedLines + 1))
	done
    echo ${RUNLIST_INFO} >> ${RUNLIST} # Print the run list info to the file
    else echo "Will not remove duplicate entries or print new entry to the file, please edit the runlist manually"
    fi
else
    echo ${RUNLIST_INFO} >> ${RUNLIST} # Print the run list info to the file
fi
