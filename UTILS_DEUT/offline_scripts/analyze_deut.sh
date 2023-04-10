#!/bin/bash

#--------------------------------------------
# Author: C. Yero
# Date: Oct 20, 2022
# email: cyero@jlab.org, cyero002@gmail.com
#--------------------------------------------

# Brief: This shell script runs the analyzer on raw replayed
#        .root files (leafs/branches) and outputs .root files
#        containing histogram objects as well as an output .txt
#        file containing a summary of the run and cuts applied. 
#        The user has the option to set/modify analysis cuts,
#        and histogram binning via an external file under inp/ dir.


# Which analysis type are we doing? "data" or "simc" (based on suffix on this shell script)
# User will need to make symlink of the original file (analyze_deut.sh) to (analyze_deut_ana_type.sh), where ana_type="simc" or "data"
ana_type=${0##*_}
ana_type=${ana_type%%.sh}  


printDataHelpMsg(){
    echo "" 
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:"
    echo ""
    echo "Usage 1 :  ./analyze_deut_${ana_type}.sh <run_number> <ana_cut> <evt_number>"
    echo ""
    echo "example: ./analyze_deut_${ana_type}.sh 16978 coin 100000"
    echo ""
    echo "------------------------------------------------------"
    echo ""
    echo "Usage 2:  ./analyze_deut_${ana_type}.sh <target> <ana_cut> <evt_number>"
    echo ""
    echo "example: ./analyze_deut_${ana_type}.sh d2 coin 100000"
    echo ""
    echo "------------------------------------------------------"
    echo ""
    echo "NOTE:"
    echo "<ana_cut> = \"bcm_calib\", \"lumi\", \"heep_singles\", \"heep_coin\", \"deep\"  "
    echo "<ana_cut> sets analysis cuts based on the type of analysis selected"
    echo ""
    echo "<target>: h2, d2, c12, dummy"
    echo ""
    echo "Usage 2: reads runs from runlist --> UTILS_DEUT/runlist/<target>_<ana_cut>.txt"
    echo ""
    echo "If no <evt_number> specified, defaults to -1 (all events) "
    
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:" 
    echo "" 
}

printSIMCHelpMsg(){
    echo "" 
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:"
    echo ""
    echo "Usage:  ./analyze_deut_${ana_type}.sh <ana_cut> <evt_number>"
    echo ""
    echo "<ana_cut> = \"heep_singles\", \"heep_coin\", \"deep\" "
    echo "<ana_cut> sets analysis cuts based on the type of analysis selected" 
    echo ""
    echo "If no <evt_number> specified, defaults to -1 (all events) " 
    echo ""
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:" 
    echo "" 
}

HCREPLAY="/work/hallc/c-deuteron/$USER/deut_offline_replay" 
echo "HCREPLAY=${HCREPLAY}"

# change to top-level directory
echo ""
echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:"
echo ""
echo "changing to the top-level directory . . ."
echo "cd ${HCREPLAY}"
cd ${HCREPLAY}
echo ""
echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:"
echo "" 

# DATA Help Message / Usage
if [ "${ana_type}" = "data" ]; then
    # check if no input are given by user
    if [ $# -eq 0 ]; then     

	# Display help output if no argument specified
	printDataHelpMsg	
	exit 0
	
	# fool-proof make sure only options: bcm_calib, lumi, optics, heep_singles, heep_coin, MF, SRC are used
    elif [ "$ana_cut" != "bcm_calib" ] || [ "$ana_cut" != "lumi" ] || [ "$ana_cut" != "heep_singles" ] ||  [ "$ana_cut" != "heep_coin" ] || [ "$ana_cut" != "deep" ]; then 
	echo "" 
    else
	printDataHelpMsg	
	exit 0
    fi
    
    # Main User Input
    run=$1     # run number
    ana_cut=$2   # Deut kinematics type, set by user:  "bcm_calib", "lumi", "heep_singles", "heep_coin",
    evt=$3     # number of events to replay (optional, but will default to all events if none specified)


    # check if 1st argument is an integer (i.e., run number, else attempt to read from runlist)
    if [ $run -eq $run ]; then
	
	
	# check if event number is specified
	if [ -z "$evt" ]; then
	    evt=-1
	    echo "No event number spedified, defaulting to evt=${evt} (all events)"
	fi
	
	# Default arguments (not required by user as input, unless the user sets as command-line arguments)
	daq_mode="coin"
	e_arm="SHMS"
	hel_flag=0
	bcm_type="BCM1"
	bcm_thrs=5.             # beam current threhsold cut > bcm_thrs [uA]
	trig_single="trig2"    # singles trigger type to apply pre-scale factor in FullWeight, i.e. hist->Scale(Ps2_factor) 
	trig_coin="trig5"      # coin. trigger type to apply pre-scale factor in FullWeight, i.e., hist->Scale(Ps5_factor)
	skim_flag=0
	
	# deut analysis script
	prod_script="UTILS_DEUT/main_analysis.cpp"
	
	# command to run analysis script ( if doing SIMC, all arguments except e_arm, ana_type and ana_cu are irrelevant )
	run_deut="root -l -q -b  \"${prod_script}( ${run},    ${evt},           
                                    \\\"${daq_mode}\\\",  \\\"${e_arm}\\\",  
                                   \\\"${ana_type}\\\", \\\"${ana_cut}\\\",
                                    ${hel_flag},                        
                                   \\\"${bcm_type}\\\", ${bcm_thrs},
                                   \\\"${trig_single}\\\", \\\"${trig_coin}\\\", ${skim_flag}
                     )\""
	{
	    echo ""
	    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
	    echo "" 
	    date
	    echo ""
	    echo ""
	    echo "Running Deut Replay Analysis on the run ${run}:"
	    echo " -> SCRIPT:  ${prod_script}"
	    echo " -> RUN:     ${run}"
	    echo " -> NEVENTS: ${evt}"
	    echo " -> COMMAND: ${run_deut}"
	    echo ""
	    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
	    
	    sleep 2
	    eval ${run_deut}
	    exit 0  
	}
	
    elif [[ -n ${run//[0-9]/} ]]; then
	# if no run-number is detected order of arguments now is as follows: target, kin, evt (run number excluded)
	# read 1st and 2nd arguments
	target=$1
	kin=$2
	evt=$3
	
	# check if event number is specified
	if [ -z "$evt" ]; then
	    evt=-1
	    echo "No event number spedified, defaulting to evt=${evt} (all events)"
	fi
	

	# runlist (when combining runs make sure they are the same kind of target and kin. type)
	filename="UTILS_DEUT/runlist/${target}_${kin}.txt"

	for run in $(cat $filename) ; do    
	    
	    
	    echo "reading run ----> ${run}" 
	    # Default arguments (not required by user as input, unless the user sets as command-line arguments)
	    daq_mode="coin"
	    e_arm="SHMS"
	    hel_flag=0
	    bcm_type="BCM1"
	    bcm_thrs=5.             # beam current threhsold cut > bcm_thrs [uA]
	    trig_single="trig2"    # singles trigger type to apply pre-scale factor in FullWeight, i.e. hist->Scale(Ps2_factor) 
	    trig_coin="trig5"      # coin. trigger type to apply pre-scale factor in FullWeight, i.e., hist->Scale(Ps5_factor)
	    skim_flag=1         # use combine runs, if a list of runs is detected
	
	    # deut analysis script
	    prod_script="UTILS_DEUT/main_analysis.cpp"
	
	    # command to run analysis script ( if doing SIMC, all arguments except e_arm, ana_type and ana_cu are irrelevant )
	    run_deut="root -l -q -b  \"${prod_script}( ${run},    ${evt},           
                                    \\\"${daq_mode}\\\",  \\\"${e_arm}\\\",  
                                   \\\"${ana_type}\\\", \\\"${ana_cut}\\\",
                                    ${hel_flag},                        
                                   \\\"${bcm_type}\\\", ${bcm_thrs},
                                   \\\"${trig_single}\\\", \\\"${trig_coin}\\\", ${skim_flag}
                     )\""

	    echo ""
	    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
	    echo "" 
	    date
	    echo ""
	    echo ""
	    echo "Running Deut Replay Analysis on the run ${run}:"
	    echo " -> SCRIPT:  ${prod_script}"
	    echo " -> RUN:     ${run}"
	    echo " -> NEVENTS: ${evt}"
	    echo " -> COMMAND: ${run_deut}"
	    echo ""
	    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
	    
	    sleep 2
	    eval ${run_deut}
	    
	    
	    
	done
    fi
fi




# SIMC Help Message / Usage
if [ "${ana_type}" = "simc" ]; then
    
    if [ -z "$2" ] || [ -z "$3" ]; then     
	
	printSIMCHelpMsg	
	exit 0
	# fool-proof make sure only options: heep_singles, heep_coin, coin are used
    elif [ [ "$ana_cut" != "heep_singles" ] ||  [ "$ana_cut" != "heep_coin" ] || [ "$ana_cut" != "deep" ] ]; then 
	echo "" 
    else
	printSIMCHelpMsg	
	exit 0
    fi

    if [ -z "$3" ]; then 
	evtNum=-1
    fi

    # Default arguments (not required by user as input, unless the user sets as command-line arguments)
    daq_mode="coin"
    e_arm="SHMS"
    hel_flag=0
    bcm_type="BCM1"
    bcm_thrs=5.             # beam current threhsold cut > bcm_thrs [uA]
    trig_single="trig2"    # singles trigger type to apply pre-scale factor in FullWeight, i.e. hist->Scale(Ps2_factor) 
    trig_coin="trig5"      # coin. trigger type to apply pre-scale factor in FullWeight, i.e., hist->Scale(Ps5_factor)
    skim_flag=0            # generate a skimmed file (file with loose coin. time cuts, centered coin. time, and 
    
    # deut analysis script
    prod_script="UTILS_DEUT/main_analysis.cpp"
    
    # command to run analysis script ( if doing SIMC, all arguments except e_arm, ana_type and ana_cu are irrelevant )
    run_deut="root -l -q -b  \"${prod_script}( ${runNum},    ${evtNum},           
            	                           \\\"${daq_mode}\\\",  \\\"${e_arm}\\\",  
                   	                    \\\"${ana_type}\\\", \\\"${ana_cut}\\\",
                            	      	        ${hel_flag},                        
                                   		\\\"${bcm_type}\\\", ${bcm_thrs},
                                   		\\\"${trig_single}\\\", \\\"${trig_coin}\\\", ${skim_flag}
                     			     )\""
    
    # Start data  analysis
    {
	
	echo "" 
	echo ""
	echo ""
	echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
	echo ""
	echo "Running Deut Data Analysis for replayed run ${runNum}:"
	echo " -> SCRIPT:  ${prod_script}"
	echo " -> RUN:     ${runNum}"
	echo " -> COMMAND: ${run_deut}"
	echo ""
	echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
	echo "" 
	echo ""
	echo ""
	eval ${run_deut}
	
    }

fi
