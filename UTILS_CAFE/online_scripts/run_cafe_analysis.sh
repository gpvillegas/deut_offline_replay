#!/bin/bash

# shell script to automatically run CaFe data replay followed by data analysis

# NOTE: During the online analysis, the user can do a partial analysis while the run
# is still going, by simply specifying the number of events, provdided that number has
# been collected already by the DAQ. For example, the user might look at 100k event replay
# in order to make count estimates, and projections on how long the run will take.


# Which analysis file type are we doing? "prod" or "sample"
ana_type=${0##*_}
ana_type=${ana_type%%.sh}


#user input
runNum=$1     # run number
kin_type=$2   # CaFe kinematics type, set by user:  "heep_singles", "heep_coin",  "MF", "SRC", depending on the production type
evtNum=$3     # number of events to replay (optional, but will default to all events if none specified)

if [ -z "$1" ] || [ -z "$2" ]; then
    echo "" 
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:"
    echo ""
    echo "Usage:  ./run_cafe_${ana_type}.sh <run_number> <kin_type> <optional evt_number>"
    echo ""
    echo "<kin_type> = \"bcm_calib\", \"lumi\", \"optics\", \"heep_singles\", \"heep_coin\", \"MF\" or \"SRC\" "
    echo ""
    echo "If you don't know which <kin_type> to choose, please ask the run coordinator ! ! ! " 
    echo ""
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:" 

    exit 0    
    # fool-proof, make sure only options: bcm_calib, lumi, optics, heep_singles, heep_coin, MF, SRC         
elif [ "$kin_type" == "bcm_calib" ] || [ "$kin_type" == "lumi" ] || [ "$kin_type" == "optics" ] || [ "$kin_type" == "heep_singles" ] ||  [ "$kin_type" == "heep_coin" ] || [ "$kin_type" == "MF" ] || [ "$kin_type" == "SRC" ]; then 
    echo ""                                                                                                                                                                                
else
    echo "" 
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:"
    echo ""
    echo "Usage: ./run_cafe_${ana_type}.sh <run_number> <kin_type> <optional evt_number>"
    echo ""     
    echo "<kin_type> = \"bcm_calib\", \"lumi\", \"optics\", \"heep_singles\", \"heep_coin\", \"MF\" or \"SRC\" "
    echo ""
    echo "If you don't know which <kin_type> to choose, please ask the run coordinator ! ! ! "   
    echo ""
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:"
    exit 0
fi

if [ -z "$3" ] && [ "${ana_type}" = "sample" ]; then
    echo "" 
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:"
    echo ""
    echo "Usage: ./run_cafe_${ana_type}.sh <run_number> <kin_type> <evt_number> <optional evt_number>"
    echo ""
    echo "No number of events was specified. Defaulting to 100k event sample"
    echo ""
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:"    
    evtNum=100000
    echo "evtNum=$evtNum"
    echo "" 
elif [ "${ana_type}" = "prod" ]; then
    echo ""
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:"
    echo ""
    echo "Usage: ./run_cafe_${ana_type}.sh <run_number> <kin_type> "
    echo ""
    echo "full event replay."
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:"
    evtNum=-1
fi

daq_mode="coin"
e_arm="SHMS"
analyze_data=1   # 1: true (analyze data), 0: false (analyze simc)
hel_flag=0
bcm_type="BCM4A"
bcm_thrs=5           # beam current threhsold cut > bcm_thrs [uA]
trig_single="trig2"    # singles trigger type to apply pre-scale factor in FullWeight, i.e. hist->Scale(Ps2_factor) 
trig_coin="trig5"      # coin. trigger type to apply pre-scale factor in FullWeight, i.e., hist->Scale(Ps5_factor)
combine_runs=0

# hcana script
if [ "${kin_type}" = "bcm_calib" ]; then
    replay_script="SCRIPTS/COIN/PRODUCTION/replay_cafe_scalers.C"
    bcm_thrs=-1      # don't apply any bcm cut 
else
    replay_script="SCRIPTS/COIN/PRODUCTION/replay_cafe.C" 
fi

# cafe serious analysis script
prod_script="UTILS_CAFE/main_analysis.cpp"
#optics_script="UTILS_CAFE/online_scripts/plotOptics.C"
optics_script="UTILS_CAFE/online_scripts/plotOptics_modified.C"   # modified by Dien Nguyen

# cafe fill run list script
fill_list_script="UTILS_CAFE/online_scripts/fill_cafe_runlist.py"

# run scripts commands
runHcana="./hcana -q \"${replay_script}(${runNum}, ${evtNum}, \\\"${ana_type}\\\")\""

runOptics="root -l -q -b \"${optics_script}(${runNum}, ${evtNum}, \\\"${ana_type}\\\")\""

runCafe="root -l -q -b \"${prod_script}( ${runNum},    ${evtNum}, 
	     	   		    \\\"${daq_mode}\\\",  \\\"${e_arm}\\\", 
				   ${analyze_data}, \\\"${kin_type}\\\", \\\"${ana_type}\\\",
          			    ${hel_flag},
                                   \\\"${bcm_type}\\\", ${bcm_thrs},
                                   \\\"${trig_single}\\\", \\\"${trig_coin}\\\", ${combine_runs}
                     )\""

fill_RunList="python ${fill_list_script} ${ana_type} ${runNum} ${evtNum}"






# Start data replay and analysis
{
    echo ""
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
    echo "" 
    date
    echo ""
    echo "Running HCANA CaFe Replay on the run ${runNum}:"
    echo " -> SCRIPT:  ${replay_script}"
    echo " -> RUN:     ${runNum}"
    echo " -> NEVENTS: ${evtNum}"
    echo " -> COMMAND: ${runHcana}"
    echo ""
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
    
    sleep 2
    eval ${runHcana}

    #--------------------------------------
    
    if [ "${kin_type}" = "optics" ]; then
	echo "" 
	echo ""
	echo ""
	echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
	echo ""
	echo "Running CaFe Optics Analysis for replayed run ${runNum}:"
	echo " -> SCRIPT:  ${optics_script}"
	echo " -> RUN:     ${runNum}"
	echo " -> COMMAND: ${runOptics}"
	echo ""
	echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
	
	sleep 2
	eval ${runOptics}
    fi
    #--------------------------------------
    echo "" 
    echo ""
    echo ""
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
    echo ""
    echo "Running CaFe Data Analysis for replayed run ${runNum}:"
    echo " -> SCRIPT:  ${prod_script}"
    echo " -> RUN:     ${runNum}"
    echo " -> COMMAND: ${runCafe}"
    echo ""
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
    
    sleep 2
    eval ${runCafe} 

    #---------------------------------------
    
    # Only full run list for production runs (i.e., full event replays)
    # sample runs (./run_cafe_sample.sh, are just for getting quick estimates to make predictions)
    if [ "${ana_type}" = "prod" ]; then
	echo "" 
	echo ""
	echo ""
	echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
	echo ""
	echo "Filling CaFe RunList for replayed run ${runNum}:"
	echo " -> SCRIPT:  ${fill_list_script}"
	echo " -> RUN:     ${runNum}"
	echo " -> COMMAND: ${fill_RunList}"
	echo ""
	echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
	
	sleep 2
	eval ${fill_RunList} 
    fi
}
