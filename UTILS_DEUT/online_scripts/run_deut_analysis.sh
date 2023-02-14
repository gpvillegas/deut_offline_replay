#!/bin/bash

# shell script to automatically run deuteron data replay followed by data analysis

# NOTE: During the online analysis, the user can do a partial analysis while the run
# is still going, by simply specifying the number of events, provdided that number has
# been collected already by the DAQ. For example, the user might look at 100k event replay
# in order to make count estimates, and projections on how long the run will take.


# Which analysis file type are we doing? "prod"
replay_type=${0##*_}
replay_type=${replay_type%%.sh}


#user input
runNum=$1     # run number
ana_cut=$2   # Deut kinematics type, set by user:  "heep_singles", "heep_coin",  "MF", "SRC", "deep",  depending on the production type
evtNum=$3     # number of events to replay (optional, but will default to all events if none specified)

if [ -z "$1" ] || [ -z "$2" ]; then
    echo "" 
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:"
    echo ""
    echo "Usage:  ./run_deut_${replay_type}.sh <run_number> <ana_cut> <evt_number>"
    echo ""
    echo "<ana_cut> = \"heep_singles\", \"heep_coin\", \"deep"\"
    echo ""
    echo "If you don't know which <ana_cut> to choose, please ask the run coordinator ! ! ! " 
    echo ""
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:" 

    exit 0    
    # fool-proof, make sure only options:  heep_singles, heep_coin, deep      
elif [ "$ana_cut" == "heep_singles" ] || [ "$ana_cut" == "heep_coin" ] || [ "$ana_cut" == "deep" ]; then 
    echo ""                                                                                                                                                                           
else
    echo "" 
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:"
    echo ""
    echo "Usage: ./run_deut_${replay_type}.sh <run_number> <ana_cut> <evt_number>"
    echo ""     
    echo "<ana_cut> = \"heep_singles\", \"heep_coin\", or \"deep\" "
    echo ""
    echo "If you don't know which <ana_cut> to choose, please ask the run coordinator ! ! ! "   
    echo ""
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:"
    exit 0
fi

if [ "${replay_type}" = "prod" ]; then
    echo ""
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:"
    echo ""
    echo "Usage: ./run_deut_${replay_type}.sh <run_number> <ana_cut> "
    echo ""
    echo "defaults to full event (-1) replay, unless event number is explicitly specified"
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:"

    if [ -z "$3" ]; then 
	evtNum=-1
    fi
fi


# check if full replay or sample events
if [[ "$evtNum" -eq -1 ]]; then
    replay_type="prod"
else
    replay_type="sample"
fi

daq_mode="coin"
e_arm="SHMS"
ana_type="data"   # "data" or "simc"
hel_flag=0
bcm_type="BCM4A"
bcm_thrs=5             # beam current threhsold cut > bcm_thrs [uA]
trig_single="trig1"    # singles trigger type to apply pre-scale factor in FullWeight, i.e. hist->Scale(Ps2_factor) 
trig_coin="trig6"      # coin. trigger type to apply pre-scale factor in FullWeight, i.e., hist->Scale(Ps5_factor)
combine_runs=0

# hcana script
if [ "${ana_cut}" = "bcm_calib" ]; then
    replay_script="SCRIPTS/COIN/PRODUCTION/replay_deut_scalers.C"
    bcm_thrs=-1      # don't apply any bcm cut 
else
    replay_script="SCRIPTS/COIN/PRODUCTION/replay_deut.C" 
fi

# deuteron serious analysis steering script
prod_script="UTILS_DEUT/main_analysis.cpp" 

# deut fill run list script
fill_list_script="UTILS_DEUT/online_scripts/fill_deut_runlist.py"  # make sure to check for deut before Feb 20/21

# run scripts commands
runHcana="./hcana -q \"${replay_script}(${runNum}, ${evtNum}, \\\"${replay_type}\\\")\""

runDeut="root -l -q -b \"${prod_script}( ${runNum},    ${evtNum}, 
	     	   		    \\\"${daq_mode}\\\",  \\\"${e_arm}\\\", 
				   \\\"${ana_type}\\\", \\\"${ana_cut}\\\",
          			    ${hel_flag},
                                   \\\"${bcm_type}\\\", ${bcm_thrs},
                                   \\\"${trig_single}\\\", \\\"${trig_coin}\\\", ${combine_runs}
                     )\""

fill_RunList="python ${fill_list_script} ${replay_type} ${runNum} ${evtNum}"






# Start data replay and analysis
{
    echo ""
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
    echo "" 
    date
    echo ""
    echo "Running HCANA d(e,e\'p) Replay on the run ${runNum}:"
    echo " -> SCRIPT:  ${replay_script}"
    echo " -> RUN:     ${runNum}"
    echo " -> NEVENTS: ${evtNum}"
    echo " -> COMMAND: ${runHcana}"
    echo ""
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
    
    sleep 2
    eval ${runHcana}

    echo "" 
    echo ""
    echo ""
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
    echo ""
    echo "Running d(e,e\'p) Data Analysis for replayed run ${runNum}:"
    echo " -> SCRIPT:  ${prod_script}"
    echo " -> RUN:     ${runNum}"
    echo " -> COMMAND: ${runDeut}"
    echo ""
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
    
    sleep 2
    eval ${runDeut} 

    #---------------------------------------
    
    # Only full run list for production runs (i.e., full event replays)
    # sample runs (./run_deut_sample.sh, are just for getting quick estimates to make predictions)
    if [ "${replay_type}" = "prod" ]; then
	echo "" 
	echo ""
	echo ""
	echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
	echo ""
	echo "Filling Deut RunList for replayed run ${runNum}:"
	echo " -> SCRIPT:  ${fill_list_script}"
	echo " -> RUN:     ${runNum}"
	echo " -> COMMAND: ${fill_RunList}"
	echo ""
	echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
	
	sleep 2
	eval ${fill_RunList} 
    fi
}
