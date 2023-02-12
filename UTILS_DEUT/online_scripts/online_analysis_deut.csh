#! /bin/tcsh

# Must be here!!
unsetenv SOURCE_SCRIPT
# Must be here!!

# C.Y. | CaFe online analysis source script

#set onlineDir="/home/cdaq/hallc-online" # This is still set to pol-he3 for now
#set onlineDir="/home/cdaq/cafe-2022"
set onlineDir="/home/cdaq/deuteron-2023"  # need to check

# Piggies love emacs
setenv EDITOR emacs

echo ""
echo "Setting up environment for the Hall-C online analysis:"

# ROOT 6 setup
#echo "--> sourcing ROOT environment"
#source $HOME/pionLT-2021/PythonPackages3.6/bin/activate.csh
#source $HOME/pionLT-2021/root-6.18.04/bin/thisroot.csh

# setenv PATH `echo $PATH | perl -pe "s[/apps/root/PRO/root/bin][]g;"`
# setenv LD_LIBRARY_PATH `echo $LD_LIBRARY_PATH | perl -pe "s[/apps/root/PRO/root/lib][]g;"`

echo "--> sourcing hcana environment"
cd "${onlineDir}/hcana"
source "${onlineDir}/hcana/setup.csh"

echo "--> sourcing replay environment"
source "${onlineDir}/deut_online_replay/setup.csh"

echo "--> changing to replay folder - ${onlineDir}/deut_online_replay"
cd "${onlineDir}/deut_online_replay"

# Add further instructions for executing main analysis script
echo ""
echo "Run './run_coin_hms.sh' to run a 50k HMS replay on the last run."
echo "Run './run_coin_shms.sh' to run a 50k SHMS replay on the last run."
echo "Run './run_deut_sample.sh' to run a 100k (default) event sample replay + physics analysis "
echo "Run './run_deut_prod.sh' to run a full event replay + physics analysis "
echo "HCANA raw data replay output is in: ./ROOTfiles and ./REPORT_OUTPUT"
echo "deuteron physics analysis output is in: ./DEUT_OUTPUT "
echo ""
