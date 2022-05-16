#!/bin/bash

# shell script to run CaFe analysis

# to pass input as command-line argument set input command you want
# to pass thru command line,  e.g., runNum=$1 target=$2,
# then this code would run as: ./run_cafe_prod $1 $2
# where $1 and $2 act as placeholders for the command line input
# the other variables are pre-defined and therefore do not require user input

runNum=$1         # run number of the input file being read in
analysis_cut=$2   # analysis cuts, set by user:  "heep",  "MF",  or "SRC", depending on the production run type

if [ -z "$1" ]; then
    echo "No run number was specified. "
    echo "e.g., ./run_cafe_sample.sh <run_number> "
   exit 0
fi

if [ -z "$2" ]; then
    echo "Specify CaFe Production Run Type: \"heep\", \"MF\" or \"SRC\"  "
    echo "e.g., ./run_cafe_sample.sh <run_number> <run_type> "
    exit 0
fi

evtNum=100000         # this is part of the input filename, 
daq_mode="coin"
e_arm="SHMS"
analysis_type="data"

hel_flag=0
bcm_type="BCM4A"
bcm_thrs=5
trig_type="trig6"
combine_runs=0

# filename base
ROOTfilePattern="OUTPUT/ROOTfiles/cafe_prod_histos_${runNum}_${evtNum}.root"

# main cafe production analysis command
CMD="root -l -q -b \"main_analysis.cpp( ${runNum},    ${evtNum}, 
	     	   		    \\\"${daq_mode}\\\",  \\\"${e_arm}\\\", 
				   \\\"${analysis_type}\\\", \\\"${analysis_cut}\\\",
          			      ${hel_flag},
                                   \\\"${bcm_type}\\\", ${bcm_thrs},
                                   \\\"${trig_type}\\\", ${combine_runs}
                     )\""


# histogram object file name
#data_ROOTfilePattern="OUTPUT/ROOTfiles/cafe_prod_histos_${runNum}_${evtNum}.root"                                                                                                                    # output_plot_fname="cafe_output_${runNum}.pdf"      

# command to make plots
#make_plots_cmd="root -l -q -b \"UTILS/make_plots.cpp(${runNum}, \\\"${data_ROOTfilePattern}\\\")\" "



echo "Running CaFe prodcution script:" $CMD
eval $CMD

#echo "Making CaFe Monitoring Plots . . ."  
#eval ${make_plots_cmd}

#echo "Opening CaFe Monitoring Plots . . ."
#openPDF="evince ${output_plot_fname}"
#eval $openPDF
