#include "baseAnalyzer.h"
#include "baseAnalyzer.cpp"
#include <iostream>

void main_simc_analysis(TString e_arm  = "SHMS", TString analysis_type = "simc", TString analysis_cut = "heep_coin")
{ 

  //----initialize baseAnalyzer (base class)----
  baseAnalyzer ba(e_arm.Data(), analysis_type.Data(), analysis_cut.Data());

  ba.run_simc_analysis();
  
}
