#! /bin/tcsh

# Must be here!!
unsetenv SOURCE_SCRIPT
# Must be here!!

# C.Y. |  NEED TO MODIFY ACCORDINGLY 

#set onlineDir="/home/cdaq/hallc-online" # This is still set to pol-he3 for now
set onlineDir="/home/cdaq/cafe-2022"

# Piggies love emacs
setenv EDITOR emacs

echo ""
echo "Setting up environment for the Hall-C online analysis:"

# ROOT 6 setup
echo "--> sourcing ROOT environment"
source $HOME/pionLT-2021/PythonPackages3.6/bin/activate.csh
source $HOME/pionLT-2021/root-6.18.04/bin/thisroot.csh

# setenv PATH `echo $PATH | perl -pe "s[/apps/root/PRO/root/bin][]g;"`
# setenv LD_LIBRARY_PATH `echo $LD_LIBRARY_PATH | perl -pe "s[/apps/root/PRO/root/lib][]g;"`

echo "--> sourcing hcana environment"
cd "${onlineDir}/hcana"
source "${onlineDir}/hcana/setup.csh"

echo "--> sourcing replay environment"
source "${onlineDir}/hallc_replay_lt/setup.csh"

echo "--> changing to replay folder - ${onlineDir}/hallc_replay_lt"
cd "${onlineDir}/hallc_replay_lt"
# Add further instructions for executing main analysis script
echo ""
echo "Run './run_coin_hms.sh' to run a 50k HMS replay on the last run."
echo "Run './run_coin_shms.sh' to run a 50k SHMS replay on the last run."
echo "Run './run_PionLT.sh' to run the physics analysis - ### Follow the instructions! ###"
echo "Outut ROOTfiles are in sub folders sorted by run type in - ROOTfiles/"
echo "Reportfiles and summaries are in sub folders sorted by run type in - REPORT_OUTPUT/"
echo ""
