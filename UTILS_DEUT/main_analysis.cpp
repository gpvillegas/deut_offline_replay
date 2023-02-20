#include "baseAnalyzer.h"
#include "baseAnalyzer.cpp"
#include <iostream>


void main_analysis(int     run           = 3243,   int evtNum           = -1,
		   TString daq_mode      = "coin", TString e_arm        = "SHMS",
		   TString analysis_type = "data", TString analysis_cut = "bcm_calib",
		   Bool_t  hel_flag     = 0, TString bcm_type  = "BCM4A",  double bcm_thrs        = 5,
		   TString trig_single = "T1", TString trig_coin = "T6",  Bool_t skim_flag    = 0		   
		   )
{


  // initialize baseAnalyzer (base class)
  baseAnalyzer ba(run, evtNum, daq_mode.Data(), e_arm.Data(), analysis_type.Data(), analysis_cut.Data(),
		  hel_flag, bcm_type.Data(), bcm_thrs, trig_single.Data(),
		  trig_coin.Data(), skim_flag);


  // data analysis
  if(analysis_type=="data"){
    
    // bcm scalers analysis
    if(analysis_cut=="bcm_calib"){
      ba.run_scalers();
    }
    // standard data analysis
    else{
      ba.run_online_data_analysis();
    }
    
  }

  // simc analysis
  if(analysis_type=="simc"){
    
    //only relevant arguments for SIMC: e_arm, analysis_type, analysis_cut
    ba.run_simc_analysis();
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
