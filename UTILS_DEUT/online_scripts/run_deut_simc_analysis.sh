#!/bin/bash

#----------------------------
# Author: C. Yero
# Sep 10, 2022
# e-mail: cyero@jlab.org,
#         cyero002@gmail.com
#----------------------------

# shell script to automatically run deuteron SIMC analysis

# NOTE: During the online analysis, the user can do SIMC analysis of
# deuteron exp. heep, deep kinematics, with assumed targets:
# hydrogen (heep) , deuteron (deep).
# The SIMC raw data file is assumed to exist, and this script will
# ONLY analyze the file (i.e, weight it, apply exact cuts as data, and


#user input
kin_type=$1   # deuteron kinematics type, set by user: "heep_coin",  "deep"


if [ -z "$1" ]; then
    echo "" 
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:"
    echo ""
    echo "Usage:  ./run_deut_simc_analysis.sh <kin_type>"
    echo ""
    echo "<kin_type> = \"heep_coin\" or \"deep\" "
    echo ""
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:" 

    exit 0    
# fool-proof, make sure only options: bcm_calib, lumi, optics, heep_singles, heep_coin, MF, SRC                                                                                         
elif [ "$kin_type" == "heep_coin" ] || [ "$kin_type" == "deep" ]; then                                                                                
    echo ""                                                                                                                                                       
else               
    echo ""                                                                                                                                                       
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:"                                                                                                         
    echo ""                                                                                                                                                                                
    echo "Usage:  ./run_deut_simc_analysis.sh <kin_type>"                                                                                                                             
    echo ""                                                                                                                                                                         
    echo "<kin_type> = \"heep_coin\" or \"deep\" "                                                                                                                  
    echo ""                                                                                                                                                                      
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:"  
    exit 0
fi

e_arm="SHMS"
analysis_type="simc"


# deuteron SIMC serious analysis script
prod_script="UTILS_DEUT/main_simc_analysis.cpp"


# run scripts commands

runDeut="root -l -q -b \"${prod_script}( \\\"${e_arm}\\\", \\\"${analysis_type}\\\", \\\"${kin_type}\\\" )\""

# Start SIMC analysis
{
  
    #--------------------------------------
    echo "" 
    echo ""
    echo ""
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
    echo ""
    echo "Running Deuteron SIMC Analysis for ${kin_type}:"
    echo " -> SCRIPT:  ${prod_script}"
    echo " -> COMMAND: ${runCafe}"
    echo ""
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
    
    sleep 2
    eval ${runDeut} 

}
