#!/bin/bash

#----------------------------
# Author: C. Yero
# Sep 10, 2022
# e-mail: cyero@jlab.org,
#         cyero002@gmail.com
#----------------------------

# shell script to automatically run CaFe SIMC analysis

# NOTE: During the online analysis, the user can do SIMC analysis of
# CaFe heep, MF or SRC kinematics, with assumed targets:
# hydrogen (heep) , carbon-12 (MF), deuterium (SRC), respectively.
# The SIMC raw data file is assumed to exist, and this script will
# ONLY analyze the file (i.e, weight it, apply exact cuts as data, and
# (in the near future), will be able to scale to other targets using target density ratio


#user input
kin_type=$1   # CaFe kinematics type, set by user: "heep_coin",  "MF" or "SRC", depending on the production type


if [ -z "$1" ]; then
    echo "" 
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:"
    echo ""
    echo "Usage:  ./run_cafe_simc_analysis.sh <kin_type>"
    echo ""
    echo "<kin_type> = \"heep_coin\", \"MF\" or \"SRC\" "
    echo ""
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:" 

    exit 0    
fi

e_arm="SHMS"
analyze_data=0   # 1: true (analyze data), 0: false (analyze simc)


# cafe SIMC serious analysis script
prod_script="UTILS_CAFE/main_simc_analysis.cpp"


# run scripts commands

runCafe="root -l -q -b \"${prod_script}( \\\"${e_arm}\\\", ${analyze_data}, \\\"${kin_type}\\\" )\""

# Start SIMC analysis
{
  
    #--------------------------------------
    echo "" 
    echo ""
    echo ""
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
    echo ""
    echo "Running CaFe SIMC Analysis for ${kin_type}:"
    echo " -> SCRIPT:  ${prod_script}"
    echo " -> COMMAND: ${runCafe}"
    echo ""
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
    
    sleep 2
    eval ${runCafe} 

}
