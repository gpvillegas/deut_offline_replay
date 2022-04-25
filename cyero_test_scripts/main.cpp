#include "baseAnalyzer.h"
#include "baseAnalyzer.cpp"
#include <iostream>

int main()
{

  /*
    What are the various inputs?
    baseAnalyzer (irun, mode, earm, type, tgt)
    
    irun: the run number
    
    mode: DAQ mode, can be either "coin" or "singles" 
    (The mode can be deremined by looking at the Start_of_Run log entry)
    
    earm: electron arm, can be either "HMS" or "SHMS"

    type: analysis type, can be either "data" or "simc"

    tgt: can be either "LH2" or "LD2" 
	  
    **NOTE: For each of these inputs, additional options can be added provided 
    they are defined as well in the baseAnalyzer.h(.cpp). See examples of the
    options above to have a better idea.

    -----------------------------------------------------------

    ***NOTE 2: Make sure to familiarize with the input files which are
   read in the analyzer code. See 'set_basic_cuts.inp' and 'set_basis_histos.inp'
    before running the code to have an idea of what is being done.

  */

  //-----------------------------------------
  //Read MAIN CONTROLS parameter input file to get initialization parameters
  //-----------------------------------------
  TString run_list, analysis_type, electron_arm,  daq_mode;
  Bool_t helicity_flag;
  
  //run_list       = trim(split(FindString("run_list", "main_controls.inp")[0], '=')[1]);  //OPTIONAL 
  analysis_type  = trim(split(FindString("analysis_type", "main_controls.inp")[0], '=')[1]);
  electron_arm   = trim(split(FindString("electron_arm", "main_controls.inp")[0], '=')[1]);
  daq_mode       = trim(split(FindString("daq_mode", "main_controls.inp")[0], '=')[1]);
  helicity_flag  = stoi(split(FindString("helicity", "main_controls.inp")[0], '=')[1]);


  //----initialize baseAnalyzer (base class)----
  baseAnalyzer ba(run, daq_mode.Data(), electron_arm.Data(), analysis_type.Data(), 0);
  ba.run_data_analysis();

  
  /*
  //------ LOOP OVER RUN LIST -----

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
    
  return 0;
}
