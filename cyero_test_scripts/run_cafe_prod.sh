#!/bin/bash

# shell script to run CaFe analysis

# to pass input as command-line argument set input command you want
# to pass thru command line,  e.g., runNum=$1 target=$2,
# then this code would run as: ./run_cafe_prod $1 $2
# where $1 and $2 act as placeholders for the command line input

runNum=$1
daq_mode="coin"
e_arm="SHMS"
analysis_type="data"
hel_flag=0
target=$2               # possible inputs: "LH2", "LD2", "Al_dummy", "B10", "B11"  
bcm_type="BCM4A"
bcm_thrs=5
trig_type="trig6"
combine_runs=0

CMD="root -l -q -b \"main_analysis.cpp(${runNum},    \\\"${daq_mode}\\\", 
                                   \\\"${e_arm}\\\", \\\"${analysis_type}\\\",
                                      ${hel_flag},   \\\"${target}\\\",
                                   \\\"${bcm_type}\\\", ${bcm_thrs},
                                   \\\"${trig_type}\\\", ${combine_runs}
                     )\""

echo $CMD
eval $CMD
