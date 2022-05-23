#!/bin/bash

# shell script to automatically run CaFe data replay followed by data analysis

# NOTE: During the online analysis, the user can do a partial analysis while the run
# is still going, by simply specifying the number of events, provdided that number has
# been collected already by the DAQ. For example, the user might look at 100k event replay
# in order to make count estimates, and projections on how long the run will take.

#user input
runNum=$1     # run number
kin_type=$2   # CaFe kinematics type, set by user:  "heep",  "MF",  or "SRC", depending on the production type
evtNum=$3     # number of events to replay (optional, but will default to all events if none specified)

if [ -z "$1" ]; then
    echo "No run number and/or run type was specified. "
    echo "e.g., ./run_cafe_sample.sh <run_number> <run_type> "
    echo "If you don't know which <run_type> to choose, please ask the run coordinator ! ! ! "
    echo "<run_type> = \"heep\", \"MF\" or \"SRC\" "
    exit 0
fi

if [ -z "$2" ]; then  
    echo "No run number and/or run type was specified. "
    echo "e.g., ./run_cafe_sample.sh <run_number> <run_type> " 
    echo "If you don't know which <run_type> to choose, please ask the run coordinator ! ! ! " 
    echo "<run_type> = \"heep\", \"MF\" or \"SRC\" " 
    exit 0  
fi

if [ -z "$3" ]; then
    echo "No number of events was specified. Defaulting to all events (-1)"
    echo "e.g., ./run_cafe_sample.sh <run_number> <run_type> <run_number>"
    evtNum=-1
fi

daq_mode="coin"
e_arm="SHMS"
analysis_type="data"
hel_flag=0
bcm_type="BCM4A"
bcm_thrs=5
trig_type="trig6"
combine_runs=0

# scripts
replay_script="SCRIPTS/COIN/PRODUCTION/replay_cafe.C"

prod_script="UTILS_CAFE/main_analysis.cpp"


# command to run scripts
runHcana="./hcana -q \"${replay_script}(${runNum}, ${evtNum}, \\\"prod\\\")\""

runCafe="root -l -q -b \"${prod_script}( ${runNum},    ${evtNum}, 
	     	   		    \\\"${daq_mode}\\\",  \\\"${e_arm}\\\", 
				   \\\"${analysis_type}\\\", \\\"${kin_type}\\\",
          			      ${hel_flag},
                                   \\\"${bcm_type}\\\", ${bcm_thrs},
                                   \\\"${trig_type}\\\", ${combine_runs}
                     )\""


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


}
