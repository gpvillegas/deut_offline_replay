#!/bin/bash

#user input
runNum=$1     # run number
evtNum=20000  # event number

# Which analysis file type are we doing? "prod" or "sample"
ana_type=${0##*_}
ana_type=${ana_type%%.sh}


replay_script="SCRIPTS/COIN/PRODUCTION/replay_cafe.C"
runHcana="./hcana -q \"${replay_script}(${runNum}, ${evtNum}, \\\"${ana_type}\\\")\""

echo $runHcana
