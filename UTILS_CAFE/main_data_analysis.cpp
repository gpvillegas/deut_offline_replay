#include "baseAnalyzer.h"
#include "baseAnalyzer.cpp"
#include <iostream>


void main_data_analysis(int     run           = 3243,   int evtNum           = -1,
		   TString daq_mode      = "coin", TString e_arm        = "SHMS",
		   Bool_t analyze_data = 1, TString analysis_cut = "bcm_calib", TString analysis_type= "prod",
		   Bool_t  hel_flag     = 0, TString bcm_type  = "BCM4A",  double bcm_thrs        = 5,
		   TString trig_single = "trig2", TString trig_coin = "trig5",  Bool_t combine_runs    = 0		   
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
    
	  
  */

  //----initialize baseAnalyzer (base class)----
  baseAnalyzer ba(run, evtNum, daq_mode.Data(), e_arm.Data(), analyze_data, analysis_cut.Data(), analysis_type.Data(), hel_flag, bcm_type.Data(), bcm_thrs, trig_single.Data(), trig_coin.Data(), combine_runs);

  if(analysis_cut=="bcm_calib"){
    ba.run_cafe_scalers();
  }
  else{
    ba.run_data_analysis();
  }


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
