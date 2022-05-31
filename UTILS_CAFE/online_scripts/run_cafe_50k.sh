#! /bin/bash


# Which spectrometer are we analyzing.
spec=${0##*_}
spec=${spec%%.sh}
SPEC=$(echo "$spec" | tr '[:lower:]' '[:upper:]')

# What is the last run number for the spectrometer.
# The pre-fix zero must be stripped because ROOT is ... well ROOT
#lastRun=$( \
#    ls raw/"${spec}"_all_*.dat raw/../raw.copiedtotape/"${spec}"_all_*.dat -R 2>/dev/null | perl -ne 'if(/0*(\d+)/) {prin#t "$1\n"}' | sort -n | tail -1 \
#)

# C.Y. Aug 26, 2021:  Commented this out since from now on, all runs that run in DAQ COIN Mode will be called shms_all_XXXXX.dat
#lastRun=$( \
#    ls raw/coin_all_*.dat raw/../raw.copiedtotape/coin_all_*.dat cache/coin_all_*.dat -R 2>/dev/null | perl -ne 'if(/0*(\d+)/) {print "$1\n"}' | sort -n | tail -1 \
#)


lastRun=$( \
	   ls raw/shms_all_*.dat raw/../raw.copiedtotape/shms_all_*.dat cache/shms_all_*.dat -R 2>/dev/null | perl -ne 'if(/0*(\d+)/) {print "$1\n"}' | sort -n | tail -1 \
       )

#
# Which run to analyze.
runNum=$1
if [ -z "$runNum" ]; then
  runNum=$lastRun
fi

# How many events to analyze.
numEvents=$2
if [ -z "$numEvents" ]; then
 numEvents=50000 
fi

# Which scripts to run.
script="SCRIPTS/COIN/PRODUCTION/replay_cafe.C"
config="CONFIG/${SPEC}/PRODUCTION/${spec}_coin_production.cfg"
expertConfig="CONFIG/${SPEC}/PRODUCTION/${spec}_coin_production_expert.cfg"

#Define some useful directories
rootFileDir="./ROOTfiles"
monRootDir="./HISTOGRAMS/${spec}50k/ROOT/"
monPdfDir="./HISTOGRAMS/${spec}50k/PDF/"
reportFileDir="./REPORT_OUTPUT/${spec}50k/"
reportMonDir="./UTIL_OL/REP_MON" 
reportMonOutDir="./MON_OUTPUT/" 

# Name of the report monitoring file
reportMonFile="reportMonitor_${spec}_${runNum}_${numEvents}.txt" 

# Which commands to run.
runHcana="./hcana -q \"${script}(${runNum}, ${numEvents}, \\\"${spec}50k\\\")\""
runOnlineGUI="./online -f ${config} -r ${runNum}"
saveOnlineGUI="./online -f ${config} -r ${runNum} -P"
saveExpertOnlineGUI="./online -f ${expertConfig} -r ${runNum} -P"
runReportMon="./${reportMonDir}/reportSummary.py ${runNum} ${numEvents} ${spec} coin"
openReportMon="emacs ${reportMonOutDir}/${reportMonFile}"  


# What is base name of onlineGUI monitoring output.                                                             
outGUI="${spec}_coin_production_${runNum}"
outGUIexpert="${spec}_coin_production_expert_${runNum}" 

outFileBase="cafe_replay_${spec}50k_monitoring_${runNum}"
outExpertFileBase="cafe_replay_${spec}50k_monitoring_expert_${runNum}"
outFileMonitor="output.txt" 

# Name of the replay ROOT file
replayFile="cafe_replay_${spec}50k_${runNum}"
rootFile="${replayFile}_${numEvents}.root"
latestRootFile="${rootFileDir}/${spec}50k/${replayFile}_latest.root"

# Names of the monitoring ROOTfile
monRootFile=${outFileBase}".root"
monExpertRootFile=${outExpertFileBase}".root" 

# Names of the monitoring PDF 
monPdfFile=${outFileBase}".pdf"
monExpertPdfFile=${outExpertFileBase}".pdf"

latestMonRootFile="cafe_replay_${spec}50k_monitoring_latest.root"
latestMonPdfFile="cafe_replay_${spec}50k_monitoring_latest.pdf" 

# Where to put log.
#reportFile="${reportFileDir}/replay_${spec}_coin_production_${runNum}_${numEvents}.txt"
reportFile="${reportFileDir}cafe_${spec}50k_${runNum}_${numEvents}.report"
#summaryFile="${reportFileDir}/summary_production_${runNum}_${numEvents}.txt"

# Replay out files
replayReport="${reportFileDir}/replayReport_${spec}_production_${runNum}_${numEvents}.txt"

# Start analysis and monitoring plots.
{
  echo ""
  echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
  echo "" 
  date
  echo ""
  echo "Running ${SPEC} COIN analysis on the run ${runNum}:"
  echo " -> SCRIPT:  ${script}"
  echo " -> RUN:     ${runNum}"
  echo " -> NEVENTS: ${numEvents}"
  echo " -> COMMAND: ${runHcana}"
  echo ""
  echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="

  sleep 2
  eval ${runHcana}

  
  # Link the ROOT file to latest for online monitoring
  ln -sf ${rootFile} ${latestRootFile}

  
  echo "" 
  echo ""
  echo ""
  echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
  echo ""
  echo "Running onlineGUI for analyzed ${SPEC} COIN run ${runNum}:"
  echo " -> CONFIG:  ${config}"
  echo " -> RUN:     ${runNum}"
  echo " -> COMMAND: ${runOnlineGUI}"
  echo ""
  echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="

  sleep 2
  cd onlineGUI
  eval ${runOnlineGUI} 
  eval ${saveOnlineGUI}
  mv "${outGUI}.pdf" "../${monPdfDir}${monPdfFile}"

  #eval ${saveExpertOnlineGUI}   # the expert onlineGUI pdf histograms will have expert-level plots + onlineGUI plots
  #mv "${outGUIexpert}.pdf" "../${monPdfDir}${monExpertPdfFile}"

  cd ../${monPdfDir}
  
  # Link onlineGUI monitoring plots to latest monitoring
  ln -sf ${monPdfFile} ${latestMonPdfFile}
  
  echo "" 
  echo ""
  echo ""
  echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
  echo ""
  echo "Done analyzing ${SPEC} run ${runNum}."
  echo ""
  echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="

  sleep 2

  echo ""
  echo ""
  echo ""
  echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
  echo ""
  echo "Generating report file monitoring data file ${SPEC} run ${runNum}."   
  echo "" 
  echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=" 


  #eval ${runReportMon}  
  #mv "${outFileMonitor}" "${reportMonOutDir}/${reportMonFile}" 
  #eval ${openReportMon}   
  #eval "emacs ${reportFile}"
  sleep 2
                                                                                        
  echo ""                                                                                                                                                                           
  echo ""                                                                                                                                                                                   
  echo ""                                                                                                                                            
  echo "-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|"                           
  echo ""                                                   
  echo "So long and thanks for all the fish!"                                             
  echo ""          
  echo "-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|"                                               
  echo ""                                                                                                                                                
  echo ""                                                                                    
  echo ""                         

} 2>&1 | tee "${replayReport}"


