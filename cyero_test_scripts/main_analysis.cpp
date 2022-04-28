#include "baseAnalyzer.h"
#include "baseAnalyzer.cpp"
#include <iostream>

void main_analysis(int     run       = 3289,       TString daq_mode      = "coin",
		   TString e_arm     = "SHMS",   TString analysis_type = "data",
		   Bool_t  hel_flag  = 0,        TString target        = "LD2",
		   TString bcm_type  = "BCM4A",  double bcm_thrs       = 5,
		   TString trig_type = "trig6",  Bool_t combine_runs   = 0
		   
		   )
{ // argumnets to add: target, bcm_type, bcm_thrs, trig_type

  /*
    What are the various inputs passed from main to baseAnalyzer?
    ----------------------
    run: the run number
    
    daq_mode: DAQ mode, can be either "coin" or "singles"  (defaults to "coin")
             (I think we'll never go back to single-arm DAQ, as singles can be taken while in coin. mode)
             (The mode can be deremined by looking at the Start_of_Run log entry)
    
    e_arm: electron arm, can be either "HMS" or "SHMS"

    analysis_type: analysis type, can be either "data" or "simc"

    hel_flag: helicity flag (0 : FALSE,  1: TRUE) to either disable or enable helicity variables readout/analysis (defaults to 0)
    
    target: liquid or solid targets: "C12" "LH2" or "LD2"  (can add more in accordance with experiment)
	  
  */

  //----initialize baseAnalyzer (base class)----
  baseAnalyzer ba(run, daq_mode.Data(), e_arm.Data(), analysis_type.Data(), hel_flag, target.Data(), bcm_type.Data(), bcm_thrs, trig_type.Data(), combine_runs);
  ba.run_data_analysis();

  /*
  //-----------------------------------------
  //Read MAIN CONTROLS parameter input file to get initialization parameters
  //-----------------------------------------
  TString run_list, analysis_type, electron_arm,  daq_mode;
  Bool_t helicity_flag;
  
  run_list       = trim(split(FindString("run_list", "main_controls.inp")[0], '=')[1]);
  analysis_type  = trim(split(FindString("analysis_type", "main_controls.inp")[0], '=')[1]);
  electron_arm   = trim(split(FindString("electron_arm", "main_controls.inp")[0], '=')[1]);
  daq_mode       = trim(split(FindString("daq_mode", "main_controls.inp")[0], '=')[1]);
  helicity_flag  = stoi(split(FindString("helicity", "main_controls.inp")[0], '=')[1]);
  */


  /*
  // ---- optional loop over multiple runs ---
  //Read list of data runs to analyze
  ifstream ifs;
  ifs.open(run_list.Data());
  string line;
  int run;

  //loop over each data run
  while(getline(ifs, line))
    {

      if (line[0]==('#')) continue;

      run = stoi(line);

      cout << "============================" << endl;
      cout << "ANALYZING RUN: " << run << endl;
      cout << "============================" << endl;

      //----initialize baseAnalyzer (base class)----
      baseAnalyzer ba(run, daq_mode.Data(), electron_arm.Data(), analysis_type.Data(), 0);
      ba.run_data_analysis();

       
    }

  */
  
}
