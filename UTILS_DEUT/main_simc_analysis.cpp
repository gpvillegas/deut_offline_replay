#include "baseAnalyzer.h"
#include "baseAnalyzer.cpp"
#include <iostream>

void main_simc_analysis(TString e_arm  = "SHMS", Bool_t analyze_data = 0, TString analysis_cut = "heep_coin")
{ 


  //----initialize baseAnalyzer (base class)----
  baseAnalyzer ba(e_arm.Data(), analyze_data, analysis_cut.Data());

  ba.run_simc_analysis();
  
}
