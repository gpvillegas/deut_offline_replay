/*
Author: Carlos Yero
email: cyero002@fiu.edu, cyero@jlab.org
Date Created: August 22, 2020
*/

#include "baseAnalyzer.h"
#include <iostream>
#include <stdio.h>
using namespace std;

//_______________________________________________________________________________
baseAnalyzer::baseAnalyzer( int irun=-1, string mode="", string earm="", string ana_type="", Bool_t hel_flag=0, string tgt_name="", string bcm_name="", double thrs=-1, string trig="", Bool_t combine_flag=0 )
  : run(irun), daq_mode(mode), e_arm_name(earm), analysis(ana_type), helicity_flag(hel_flag), tgt_type(tgt_name), bcm_type(bcm_name), bcm_thrs(thrs), trig_type(trig), combine_runs_flag(combine_flag)   //initialize member list 
{

  cout << "Calling BaseConstructor " << endl;
  
  //Set prefix depending on DAQ mode and electron arm (used for naming leaf variables)
  if(daq_mode=="coin" && e_arm_name=="SHMS"){
    eArm = "P";
    hArm = "H";
    e_arm = "p";
    h_arm = "h";
    nroc = "2";
    daq = "coin";    //string is used for storing either "shms", "hms" or "coin" in trigger leaf variable (T.*)  
    scl_tree_name = "TSP";  //By default, if in coin mode, Read SHMS Scaler Tree
    h_arm_name = "HMS";
  }
  
  else if(daq_mode=="coin" && e_arm_name=="HMS"){
    eArm = "H";
    hArm = "P";
    e_arm = "h";
    h_arm = "p";
    nroc = "2";
    daq = "coin";    //used for storing either "shms", "hms" or "coin" in trigger leaf variable (T.*)
    scl_tree_name = "TSP";  //By default, if in coin mode, Read SHMS Scaler Tree
    h_arm_name = "SHMS";
  }

  else if(daq_mode=="singles"){
    if(e_arm_name=="HMS") { eArm = "H"; e_arm = "h"; nroc = "1"; daq = "hms";  scl_tree_name = "TSH";}
    else {eArm = "P"; e_arm = "p", nroc = "2", daq = "shms"; scl_tree_name = "TSP";}
  }

  

  //Initialize TFile Pointers
  inROOT  = NULL;
  outROOT = NULL;
  
  //Initialize TTree Pointers
  tree        = NULL;
  scaler_tree = NULL;

  //Initialize Scaler Pointer
  evt_flag_bcm = NULL;   //scaler read array of 0 or 1 (determine if scaler read passed bcm cut)
  scal_evt_num = NULL;

  //Initialize TList Pointers
  pid_HList = NULL;
  kin_HList = NULL;
  accp_HList = NULL;
  
  //-----Initialize Histogram Pointers-----

  // Dummy histograms to store the ith and cumulative histograms (See CombineHistos() Method)
  //1D
  h_total = NULL;         //dummy histo to store hsitogram sum
  h_i     = NULL;         //dummy histo to store ith histogram from list
  //2D
  h2_total = NULL;        //dummy histo to store hsitogram sum
  h2_i     = NULL;       //dummy histo to store ith histogram from list

  
  //---------------------------------------------------------------
  //Detector Histograms (DATA ONLY) --- PID / TRACKING EFFICIENCY
  //---------------------------------------------------------------
  
  //Coincidence Time
  H_ep_ctime  = NULL;
  H_eK_ctime  = NULL;
  H_ePi_ctime = NULL;
  
  //-HMS-
  H_hCerNpeSum      = NULL;  
  H_hCalEtotNorm    = NULL;
  H_hCalEtotTrkNorm = NULL;
  H_hHodBetaNtrk    = NULL;   
  H_hHodBetaTrk    = NULL;
		  
  //-SHMS-
  H_pNGCerNpeSum    = NULL;
  H_pHGCerNpeSum    = NULL;
  H_pAeroNpeSum     = NULL;
  H_pCalEtotNorm    = NULL;
  H_pCalEtotTrkNorm = NULL;
  H_pHodBetaNtrk    = NULL;   
  H_pHodBetaTrk     = NULL;
  
  //HMS 2D PID
  H_hcal_vs_hcer     = NULL;  
  
  //sHMS 2D PID
  H_pcal_vs_phgcer   = NULL;  
  H_pcal_vs_pngcer   = NULL;  
  H_pcal_vs_paero    = NULL;   
  H_paero_vs_phgcer  = NULL; 
  H_paero_vs_pngcer  = NULL; 
  H_pngcer_vs_phgcer = NULL;

  //DATA/SIMC Histograms (MUST BE THE EXACT SAME HSITOGRAM BINNING)

  //-------------------------
  //  Kinematics Histograms
  //-------------------------

  //Primary (electron) Kinematics
  H_the      = NULL; 
  H_W        = NULL; 
  H_W2       = NULL; 
  H_Q2       = NULL; 
  H_xbj      = NULL; 
  H_nu       = NULL; 
  H_q        = NULL; 
  H_qx       = NULL; 
  H_qy       = NULL; 
  H_qz       = NULL; 
  H_thq      = NULL;  
  H_phq      = NULL;   
  H_epsilon  = NULL; 
  
  //Secondary (Hadron) Kinematics
  H_Em         = NULL;	
  H_Em_nuc     = NULL;	
  H_Pm         = NULL;	
  H_Pmx_lab    = NULL;	
  H_Pmy_lab    = NULL;	
  H_Pmz_lab    = NULL;	
  H_Pmx_q      = NULL;	
  H_Pmy_q      = NULL;	
  H_Pmz_q      = NULL;	
  H_Tx         = NULL;	
  H_Tr         = NULL;	
  H_MM         = NULL;	
  H_thxq       = NULL;	
  H_thrq       = NULL;	
  H_phxq       = NULL;	
  H_phrq       = NULL;	
  H_Tx_cm      = NULL;	
  H_Tr_cm      = NULL;	
  H_thxq_cm    = NULL;	
  H_thrq_cm    = NULL;	
  H_phxq_cm    = NULL;	
  H_phrq_cm    = NULL;	
  H_Ttot_cm    = NULL;	
  H_MandelS    = NULL;	
  H_MandelT    = NULL;	
  H_MandelU    = NULL;  

  //Kinematics Defined in HCANA (which are not in primary/secondary modules)
  H_kf    = NULL;
  H_Pf    = NULL;

  //Additional Kinematics (User-defined)
  H_thx    = NULL;
  H_MM2    = NULL;

  //----(Cosine, Sine) Histos of detecte AND recoil angles-----
  //(No need to define bins, as they are limited to -1 to 1). Will explicitly define bine in code
  //cth: cos(theta),  sth: sin(theta), cphi: cos(phi), sphi: sin(phi)
  H_cth_xq      = NULL;	 
  H_cth_rq      = NULL;	 
  H_sth_xq      = NULL;	 
  H_sth_rq      = NULL;	 
  H_cphi_xq     = NULL;	 
  H_cphi_rq     = NULL;	 
  H_sphi_xq     = NULL;	 
  H_sphi_rq     = NULL;	 
			 
  H_cth_xq_cm   = NULL;	 
  H_cth_rq_cm   = NULL;	 
  H_sth_xq_cm   = NULL;	 
  H_sth_rq_cm   = NULL;	 
  H_cphi_xq_cm  = NULL;	 
  H_cphi_rq_cm  = NULL;	 
  H_sphi_xq_cm  = NULL;	 
  H_sphi_rq_cm  = NULL;  
  
  //-------------------------
  //  Acceptance Histograms
  //-------------------------
  
   //Electron Arm Focal Plane Quantities 
  H_exfp     = NULL;   
  H_expfp    = NULL;  
  H_eyfp     = NULL;   
  H_eypfp    = NULL;  

  //Electron Arm Reconstructed Quantities
  H_eytar     = NULL;  
  H_eyptar    = NULL; 
  H_exptar    = NULL; 
  H_edelta    = NULL; 

  //Hadron Arm Focal Plane / Reconstructed Quantities
  H_hxfp     = NULL;   
  H_hxpfp    = NULL;  
  H_hyfp     = NULL;   
  H_hypfp    = NULL;  
  		      
  H_hytar     = NULL;  
  H_hyptar    = NULL; 
  H_hxptar    = NULL; 
  H_hdelta    = NULL; 
  
  //Target Quantities (tarx, tary, tarz) in Hall Coord. System 
  H_htar_x    = NULL;  
  H_htar_y    = NULL;  
  H_htar_z    = NULL;  
  		      
  H_etar_x    = NULL;  
  H_etar_y    = NULL;  
  H_etar_z    = NULL;   

  //Additional Acceptance Quantities (User-defined)
  H_ztar_diff    = NULL; 

  //Collimator Quantities (in spectrometer corrdinate system)
  H_hXColl    = NULL;  
  H_hYColl    = NULL;  
  H_eXColl    = NULL;  
  H_eYColl    = NULL;  
  
  //--2D Acceptance Histos--

  //Collimator Shape
  H_hXColl_vs_hYColl    = NULL;
  H_eXColl_vs_eYColl    = NULL;
  //Hour-Glass Shape
  H_hxfp_vs_hyfp    = NULL;
  H_exfp_vs_eyfp    = NULL;
  
}

//_______________________________________________________________________________
baseAnalyzer::~baseAnalyzer()
{
  cout << "Calling BaseDestructor " << endl;

  //Delete FileName Pointers
  delete inROOT; inROOT   = NULL;
  delete outROOT; outROOT = NULL;

  //Delete Scaler related event flag
  delete [] evt_flag_bcm; evt_flag_bcm = NULL;
  delete [] scal_evt_num; scal_evt_num = NULL;

  //Delete TList Pointers
  delete pid_HList; pid_HList = NULL;
  delete kin_HList; kin_HList = NULL;
  delete accp_HList; accp_HList = NULL;

  
  //-----------------------------
  // Detector Histogram Pointers
  //-----------------------------

  //-Coin. Time-
  delete H_ep_ctime; H_ep_ctime   = NULL;
  delete H_eK_ctime; H_eK_ctime   = NULL;
  delete H_ePi_ctime; H_ePi_ctime = NULL;
  
  //-HMS-
  delete H_hCerNpeSum;      H_hCerNpeSum      = NULL;  
  delete H_hCalEtotNorm;    H_hCalEtotNorm    = NULL;
  delete H_hCalEtotTrkNorm; H_hCalEtotTrkNorm = NULL;
  delete H_hHodBetaNtrk;    H_hHodBetaNtrk    = NULL;   
  delete H_hHodBetaTrk;     H_hHodBetaTrk     = NULL;
		  
  //-SHMS-
  delete H_pNGCerNpeSum;    H_pNGCerNpeSum    = NULL;
  delete H_pHGCerNpeSum;    H_pHGCerNpeSum    = NULL;
  delete H_pAeroNpeSum;     H_pAeroNpeSum     = NULL; 
  delete H_pCalEtotNorm;    H_pCalEtotNorm    = NULL;
  delete H_pCalEtotTrkNorm; H_pCalEtotTrkNorm = NULL;
  delete H_pHodBetaNtrk;    H_pHodBetaNtrk    = NULL;   
  delete H_pHodBetaTrk;     H_pHodBetaTrk     = NULL;   
  
  //HMS 2D PID                      
  delete H_hcal_vs_hcer;         H_hcal_vs_hcer     = NULL;
  		     	     			    
  //sHMS 2D PID	     	     		    
  delete H_pcal_vs_phgcer;   	 H_pcal_vs_phgcer   = NULL;
  delete H_pcal_vs_pngcer;   	 H_pcal_vs_pngcer   = NULL;
  delete H_pcal_vs_paero;    	 H_pcal_vs_paero    = NULL;
  delete H_paero_vs_phgcer;  	 H_paero_vs_phgcer  = NULL;
  delete H_paero_vs_pngcer;  	 H_paero_vs_pngcer  = NULL;
  delete H_pngcer_vs_phgcer; 	 H_pngcer_vs_phgcer = NULL;
  
  //Delete DATA/SIMC Histogram Pointers 

  //----------------------------------
  //   Kinematics Histograms Pointers
  //----------------------------------
  
  //Primary (electron) Kinematics
  delete H_the;          H_the      = NULL; 
  delete H_W;	         H_W        = NULL; 
  delete H_W2;	         H_W2       = NULL; 
  delete H_Q2;	         H_Q2       = NULL; 
  delete H_xbj;	         H_xbj      = NULL; 
  delete H_nu;	         H_nu       = NULL; 
  delete H_q;	         H_q        = NULL; 
  delete H_qx;	         H_qx       = NULL; 
  delete H_qy;	         H_qy       = NULL; 
  delete H_qz;	         H_qz       = NULL; 
  delete H_thq;          H_thq      = NULL; 
  delete H_phq;          H_phq      = NULL; 
  delete H_epsilon;      H_epsilon  = NULL; 
  
  //Secondary (Hadron) Kinematics
  delete H_Em;             H_Em         = NULL;	
  delete H_Em_nuc;	   H_Em_nuc     = NULL;	
  delete H_Pm;		   H_Pm         = NULL;	
  delete H_Pmx_lab;	   H_Pmx_lab    = NULL;	
  delete H_Pmy_lab;	   H_Pmy_lab    = NULL;	
  delete H_Pmz_lab;	   H_Pmz_lab    = NULL;	
  delete H_Pmx_q;	   H_Pmx_q      = NULL;	
  delete H_Pmy_q;	   H_Pmy_q      = NULL;	
  delete H_Pmz_q;	   H_Pmz_q      = NULL;	
  delete H_Tx;		   H_Tx         = NULL;	
  delete H_Tr;		   H_Tr         = NULL;	
  delete H_MM;		   H_MM         = NULL;	
  delete H_thxq;	   H_thxq       = NULL;	
  delete H_thrq;	   H_thrq       = NULL;	
  delete H_phxq;	   H_phxq       = NULL;	
  delete H_phrq;	   H_phrq       = NULL;	
  delete H_Tx_cm;	   H_Tx_cm      = NULL;	
  delete H_Tr_cm;	   H_Tr_cm      = NULL;	
  delete H_thxq_cm;	   H_thxq_cm    = NULL;	
  delete H_thrq_cm;	   H_thrq_cm    = NULL;	
  delete H_phxq_cm;	   H_phxq_cm    = NULL;	
  delete H_phrq_cm;	   H_phrq_cm    = NULL;	
  delete H_Ttot_cm;	   H_Ttot_cm    = NULL;	
  delete H_MandelS;	   H_MandelS    = NULL;	
  delete H_MandelT;	   H_MandelT    = NULL;	
  delete H_MandelU;	   H_MandelU    = NULL;  

  //Kinematics Defined in HCANA (which are not in primary/secondary modules)
  delete H_kf;   H_kf    = NULL;
  delete H_Pf;   H_Pf    = NULL;

  //Additional Kinematics (User-defined)
  delete H_thx;   H_thx    = NULL;
  delete H_MM2;   H_MM2    = NULL;

  //----(Cosine, Sine) Histos of detecte AND recoil angles-----
  //(No need to define bins, as they are limited to -1 to 1). Will explicitly define bine in code
  //cth: cos(theta),  sth: sin(theta), cphi: cos(phi), sphi: sin(phi)
  delete H_cth_xq;          H_cth_xq      = NULL;	 
  delete H_cth_rq;	    H_cth_rq      = NULL;	 
  delete H_sth_xq;	    H_sth_xq      = NULL;	 
  delete H_sth_rq;	    H_sth_rq      = NULL;	 
  delete H_cphi_xq;	    H_cphi_xq     = NULL;	 
  delete H_cphi_rq;	    H_cphi_rq     = NULL;	 
  delete H_sphi_xq;	    H_sphi_xq     = NULL;	 
  delete H_sphi_rq;	    H_sphi_rq     = NULL;	 
			  			 
  delete H_cth_xq_cm;	    H_cth_xq_cm   = NULL;	 
  delete H_cth_rq_cm;	    H_cth_rq_cm   = NULL;	 
  delete H_sth_xq_cm;	    H_sth_xq_cm   = NULL;	 
  delete H_sth_rq_cm;	    H_sth_rq_cm   = NULL;	 
  delete H_cphi_xq_cm;	    H_cphi_xq_cm  = NULL;	 
  delete H_cphi_rq_cm;	    H_cphi_rq_cm  = NULL;	 
  delete H_sphi_xq_cm;	    H_sphi_xq_cm  = NULL;	 
  delete H_sphi_rq_cm;	    H_sphi_rq_cm  = NULL;  
  
  //------------------------------------------------
  
  //-----------------------------
  //   Acceptance Histograms
  //-----------------------------
  
  //Electron Arm Focal Plane Quantities 
  delete H_exfp;       H_exfp     = NULL;  
  delete H_expfp;      H_expfp    = NULL;  
  delete H_eyfp;       H_eyfp     = NULL;  
  delete H_eypfp;      H_eypfp    = NULL;  

  //Electron Arm Reconstructed Quantities
  delete H_eytar;      H_eytar     = NULL; 
  delete H_eyptar;     H_eyptar    = NULL; 
  delete H_exptar;     H_exptar    = NULL; 
  delete H_edelta;     H_edelta    = NULL; 

  //Hadron Arm Focal Plane / Reconstructed Quantities
  delete H_hxfp;       H_hxfp     = NULL;  
  delete H_hxpfp;      H_hxpfp    = NULL;  
  delete H_hyfp;       H_hyfp     = NULL;  
  delete H_hypfp;      H_hypfp    = NULL;  
  		        		      
  delete H_hytar;      H_hytar     = NULL; 
  delete H_hyptar;     H_hyptar    = NULL; 
  delete H_hxptar;     H_hxptar    = NULL; 
  delete H_hdelta;     H_hdelta    = NULL; 
  
  //Target Quantities (tarx, tary, tarz) in Hall Coord. System 
  delete H_htar_x;      H_htar_x    = NULL; 
  delete H_htar_y;      H_htar_y    = NULL; 
  delete H_htar_z;      H_htar_z    = NULL; 
  		        		      
  delete H_etar_x;      H_etar_x    = NULL; 
  delete H_etar_y;      H_etar_y    = NULL; 
  delete H_etar_z;      H_etar_z    = NULL; 

  //Additional Acceptance Quantities (User-defined)
  delete H_ztar_diff;   H_ztar_diff    = NULL; 
  
  //Collimator Quantities (in spectrometer corrdinate system)
  delete H_hXColl;     H_hXColl    = NULL;  
  delete H_hYColl;     H_hYColl    = NULL;  
  delete H_eXColl;     H_eXColl    = NULL;  
  delete H_eYColl;     H_eYColl    = NULL;  
  
  //--2D Acceptance Histos--

  //Collimator Shape
  delete H_hXColl_vs_hYColl;    H_hXColl_vs_hYColl    = NULL;
  delete H_eXColl_vs_eYColl;    H_eXColl_vs_eYColl    = NULL;
  //Hour-Glass Shape
  delete H_hxfp_vs_hyfp;        H_hxfp_vs_hyfp    = NULL;
  delete H_exfp_vs_eyfp;        H_exfp_vs_eyfp    = NULL;

  
}

//_______________________________________________________________________________
void baseAnalyzer::ReadInputFile()
{
  cout << "Calling Base ReadInputFiles() . . . " << endl;
  
  /*
    Brief: Read Input Files: either main control parameters, tracking eff. cuts
    or data analysis cuts. This method allows portability and reusability of the 
    parameters to be read, as well as makes parameter file modifications easier
    for any class that inherits fromt this class. The user inherits this method
    and can easily modify or even overwrite the format of how parameters are read,
    in accordance with the analysis requirements.
   */

  input_FileNamePattern = "inp/set_basic_filenames.inp";
  input_CutFileName     = "inp/set_basic_cuts.inp";
  input_HBinFileName    = "inp/set_basic_histos.inp";

  //==========================================
  //     READ FILE NAME PATTERN
  //==========================================

  TString temp; //temporary string placeholder


  
  //-------------------------------
  //----INPUTS (USER READS IN)-----
  //-------------------------------
  
  //Define Input (.root) File Name Patterns (read principal ROOTfile from experiment)
  temp = trim(split(FindString("input_ROOTfilePattern", input_FileNamePattern.Data())[0], '=')[1]);
  data_InputFileName = Form(temp.Data(),  run);


  //Define Input (.report) File Name Pattern (read principal REPORTfile from experiment)
  temp = trim(split(FindString("input_REPORTPattern", input_FileNamePattern.Data())[0], '=')[1]);
  data_InputReport = Form(temp.Data(), run);


  //----------------------------------
  //----OUTPUTS (USER WRITES OUT)-----
  //----------------------------------

  //Define Output (.root) File Name Pattern (analyzed histos are written to this file)
  temp = trim(split(FindString("output_ROOTfilePattern", input_FileNamePattern.Data())[0], '=')[1]);
  data_OutputFileName = Form(temp.Data(), run);

  //Define Output (.root) File Name Pattern (analyzed combined histos are written to this file)
  temp = trim(split(FindString("output_ROOTfilePattern_final", input_FileNamePattern.Data())[0], '=')[1]);
  data_OutputFileName_combined = temp.Data();

  //Define Output (.txt) File Name Pattern (analysis report is written to this file)
  temp = trim(split(FindString("output_REPORTPattern", input_FileNamePattern.Data())[0], '=')[1]);
  report_OutputFileName = temp.Data();


  //==========================================
  //     READ TRACKING EFFICIENCY CUTS
  //==========================================
  
  //HMS Tracking Efficiency Cut Flags / Limits
  hdc_ntrk_cut_flag = stoi(split(FindString("hdc_ntrk_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_hdc_ntrk_min = stod(split(FindString("c_hdc_ntrk_min", input_CutFileName.Data())[0], '=')[1]);
  
  hScinGood_cut_flag = stoi(split(FindString("hScinGood_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  
  hcer_cut_flag = stoi(split(FindString("hcer_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_hnpeSum_min = stod(split(FindString("c_hnpeSum_min", input_CutFileName.Data())[0], '=')[1]);
  c_hnpeSum_max = stod(split(FindString("c_hnpeSum_max", input_CutFileName.Data())[0], '=')[1]);
  
  hetotnorm_cut_flag = stoi(split(FindString("hetotnorm_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_hetotnorm_min = stod(split(FindString("c_hetotnorm_min", input_CutFileName.Data())[0], '=')[1]);
  c_hetotnorm_max = stod(split(FindString("c_hetotnorm_max", input_CutFileName.Data())[0], '=')[1]);
  
  hBeta_notrk_cut_flag = stoi(split(FindString("hBeta_notrk_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_hBetaNtrk_min = stod(split(FindString("c_hBetaNtrk_min", input_CutFileName.Data())[0], '=')[1]);
  c_hBetaNtrk_max = stod(split(FindString("c_hBetaNtrk_max", input_CutFileName.Data())[0], '=')[1]);
  
  //SHMS Tracking Efficiency Cut Flags / Limits
  pdc_ntrk_cut_flag = stoi(split(FindString("pdc_ntrk_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_pdc_ntrk_min = stod(split(FindString("c_pdc_ntrk_min", input_CutFileName.Data())[0], '=')[1]);
  
  pScinGood_cut_flag = stoi(split(FindString("pScinGood_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  
  pngcer_cut_flag = stoi(split(FindString("pngcer_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_pngcer_npeSum_min = stod(split(FindString("c_pngcer_npeSum_min", input_CutFileName.Data())[0], '=')[1]);
  c_pngcer_npeSum_max = stod(split(FindString("c_pngcer_npeSum_max", input_CutFileName.Data())[0], '=')[1]);

  phgcer_cut_flag = stoi(split(FindString("phgcer_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_phgcer_npeSum_min = stod(split(FindString("c_phgcer_npeSum_min", input_CutFileName.Data())[0], '=')[1]);
  c_phgcer_npeSum_max = stod(split(FindString("c_phgcer_npeSum_max", input_CutFileName.Data())[0], '=')[1]);
  
  petotnorm_cut_flag = stoi(split(FindString("petotnorm_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_petotnorm_min = stod(split(FindString("c_petotnorm_min", input_CutFileName.Data())[0], '=')[1]);
  c_petotnorm_max = stod(split(FindString("c_petotnorm_max", input_CutFileName.Data())[0], '=')[1]);
  
  pBeta_notrk_cut_flag = stoi(split(FindString("pBeta_notrk_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_pBetaNtrk_min = stod(split(FindString("c_pBetaNtrk_min", input_CutFileName.Data())[0], '=')[1]);
  c_pBetaNtrk_max = stod(split(FindString("c_pBetaNtrk_max", input_CutFileName.Data())[0], '=')[1]);
  
  
  //==========================================
  
  
  
  
  //==========================================
  //     READ Data/SIMC ANALYSIS CUTS
  //==========================================
  
  //------PID Cuts-----
  
  //Coincidence time cuts (check which coin. time cut is actually being applied. By default: electron-proton cut is being applied)
  
  eKctime_pidCut_flag = stoi(split(FindString("eKctime_pidCut_flag", input_CutFileName.Data())[0], '=')[1]);  
  cpid_eKctime_min = stod(split(FindString("cpid_eKctime_min", input_CutFileName.Data())[0], '=')[1]);
  cpid_eKctime_max = stod(split(FindString("cpid_eKctime_max", input_CutFileName.Data())[0], '=')[1]);
  
  ePictime_pidCut_flag = stoi(split(FindString("ePictime_pidCut_flag", input_CutFileName.Data())[0], '=')[1]);
  cpid_ePictime_min = stod(split(FindString("cpid_ePictime_min", input_CutFileName.Data())[0], '=')[1]);
  cpid_ePictime_max = stod(split(FindString("cpid_ePictime_max", input_CutFileName.Data())[0], '=')[1]);
  
  ePctime_pidCut_flag = stoi(split(FindString("ePctime_pidCut_flag", input_CutFileName.Data())[0], '=')[1]);
  cpid_ePctime_min = stod(split(FindString("cpid_ePctime_min", input_CutFileName.Data())[0], '=')[1]);
  cpid_ePctime_max = stod(split(FindString("cpid_ePctime_max", input_CutFileName.Data())[0], '=')[1]);
  
  //(SHMS PID) Calorimeter Total Energy Normalized By Track Momentum
  petot_trkNorm_pidCut_flag = stoi(split(FindString("petot_trkNorm_pidCut_flag", input_CutFileName.Data())[0], '=')[1]);
  cpid_petot_trkNorm_min = stod(split(FindString("cpid_petot_trkNorm_min", input_CutFileName.Data())[0], '=')[1]);
  cpid_petot_trkNorm_max = stod(split(FindString("cpid_petot_trkNorm_max", input_CutFileName.Data())[0], '=')[1]);
  
  //(SHMS PID) Noble Gas Cherenkov
  pngcer_pidCut_flag = stoi(split(FindString("pngcer_pidCut_flag", input_CutFileName.Data())[0], '=')[1]);
  cpid_pngcer_npeSum_min = stod(split(FindString("cpid_pngcer_npeSum_min", input_CutFileName.Data())[0], '=')[1]);
  cpid_pngcer_npeSum_max = stod(split(FindString("cpid_pngcer_npeSum_max", input_CutFileName.Data())[0], '=')[1]);
  
  //(SHMS PID) Heavy Gas Cherenkov
  phgcer_pidCut_flag = stoi(split(FindString("phgcer_pidCut_flag", input_CutFileName.Data())[0], '=')[1]);
  cpid_phgcer_npeSum_min = stod(split(FindString("cpid_phgcer_npeSum_min", input_CutFileName.Data())[0], '=')[1]);
  cpid_phgcer_npeSum_max = stod(split(FindString("cpid_phgcer_npeSum_max", input_CutFileName.Data())[0], '=')[1]);
  
  //(SHMS PID) Aerogel Cherenkov
  paero_pidCut_flag = stoi(split(FindString("paero_pidCut_flag", input_CutFileName.Data())[0], '=')[1]);
  cpid_paero_npeSum_min = stod(split(FindString("cpid_paero_npeSum_min", input_CutFileName.Data())[0], '=')[1]);
  cpid_paero_npeSum_max = stod(split(FindString("cpid_paero_npeSum_max", input_CutFileName.Data())[0], '=')[1]);
  
  //(HMS PID) Calorimeter Total Energy Normalized By Track Momentum
  hetot_trkNorm_pidCut_flag = stoi(split(FindString("hetot_trkNorm_pidCut_flag", input_CutFileName.Data())[0], '=')[1]);
  cpid_hetot_trkNorm_min = stod(split(FindString("cpid_hetot_trkNorm_min", input_CutFileName.Data())[0], '=')[1]);
  cpid_hetot_trkNorm_max = stod(split(FindString("cpid_hetot_trkNorm_max", input_CutFileName.Data())[0], '=')[1]);
  
  //(HMS PID) Gas Cherenkov
  hcer_pidCut_flag = stoi(split(FindString("hcer_pidCut_flag", input_CutFileName.Data())[0], '=')[1]);
  cpid_hcer_npeSum_min = stod(split(FindString("cpid_hcer_npeSum_min", input_CutFileName.Data())[0], '=')[1]);
  cpid_hcer_npeSum_max = stod(split(FindString("cpid_hcer_npeSum_max", input_CutFileName.Data())[0], '=')[1]);
  
  //-----Kinematics Cuts------
  
  //4-Momentum Transfers [GeV^2]
  Q2_cut_flag = stoi(split(FindString("Q2_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_Q2_min = stod(split(FindString("c_Q2_min", input_CutFileName.Data())[0], '=')[1]);
  c_Q2_max = stod(split(FindString("c_Q2_max", input_CutFileName.Data())[0], '=')[1]);
  
  //Missing Energy [GeV]
  Em_cut_flag = stoi(split(FindString("Em_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_Em_min = stod(split(FindString("c_Em_min", input_CutFileName.Data())[0], '=')[1]);
  c_Em_max = stod(split(FindString("c_Em_max", input_CutFileName.Data())[0], '=')[1]);
  
  //Invariant Mass, W [GeV]
  W_cut_flag = stoi(split(FindString("W_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_W_min = stod(split(FindString("c_W_min", input_CutFileName.Data())[0], '=')[1]);
  c_W_max = stod(split(FindString("c_W_max", input_CutFileName.Data())[0], '=')[1]);
  
  //Missing Mass Cut (Check which MM Cut is actually being applied: By default, it should be proton MM)
  
  //Kaons
  MM_K_cut_flag = stoi(split(FindString("MM_K_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_MM_K_min = stod(split(FindString("c_MM_K_min", input_CutFileName.Data())[0], '=')[1]);
  c_MM_K_max = stod(split(FindString("c_MM_K_max", input_CutFileName.Data())[0], '=')[1]);
  
  //Pions
  MM_Pi_cut_flag = stoi(split(FindString("MM_Pi_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_MM_Pi_min = stod(split(FindString("c_MM_Pi_min", input_CutFileName.Data())[0], '=')[1]);
  c_MM_Pi_max = stod(split(FindString("c_MM_Pi_max", input_CutFileName.Data())[0], '=')[1]);
  
  //Protons
  MM_P_cut_flag = stoi(split(FindString("MM_P_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_MM_P_min = stod(split(FindString("c_MM_P_min", input_CutFileName.Data())[0], '=')[1]);
  c_MM_P_max = stod(split(FindString("c_MM_P_max", input_CutFileName.Data())[0], '=')[1]);
  
  //------Acceptance Cuts-------
  
  //Hadron Arm
  hdelta_cut_flag = stoi(split(FindString("hdelta_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_hdelta_min = stod(split(FindString("c_hdelta_min", input_CutFileName.Data())[0], '=')[1]);
  c_hdelta_max = stod(split(FindString("c_hdelta_max", input_CutFileName.Data())[0], '=')[1]);
  
  //Electron Arm
  edelta_cut_flag = stoi(split(FindString("edelta_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_edelta_min = stod(split(FindString("c_edelta_min", input_CutFileName.Data())[0], '=')[1]);
  c_edelta_max = stod(split(FindString("c_edelta_max", input_CutFileName.Data())[0], '=')[1]);
  
  // Z-Reaction Vertex Difference Cut
  ztarDiff_cut_flag = stoi(split(FindString("ztarDiff_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_ztarDiff_min = stod(split(FindString("c_ztarDiff_min", input_CutFileName.Data())[0], '=')[1]);
  c_ztarDiff_max = stod(split(FindString("c_ztarDiff_max", input_CutFileName.Data())[0], '=')[1]);
  
 
  
}
  

//_______________________________________________________________________________
void baseAnalyzer::ReadReport()
{

  //Brief: Read Necessary Quantities from Report File
    
  cout << "Calling Base ReadReport() " << endl;
  
  string temp;

  //Read Pre-Scale Factors
  temp =  FindString("Ps1_factor", data_InputReport)[0]; 
  Ps1_factor = stod(split(temp, '=')[1]);

  temp =  FindString("Ps2_factor", data_InputReport)[0]; 
  Ps2_factor = stod(split(temp, '=')[1]);
  
  temp =  FindString("Ps3_factor", data_InputReport)[0]; 
  Ps3_factor = stod(split(temp, '=')[1]);
  
  temp =  FindString("Ps4_factor", data_InputReport)[0]; 
  Ps4_factor = stod(split(temp, '=')[1]);
  
  temp =  FindString("Ps5_factor", data_InputReport)[0]; 
  Ps5_factor = stod(split(temp, '=')[1]);
  
  temp =  FindString("Ps6_factor", data_InputReport)[0]; 
  Ps6_factor = stod(split(temp, '=')[1]);
  

  
}

//_______________________________________________________________________________
void baseAnalyzer::SetHistBins()
{
  cout << "Calling Base SetHistBins()  " << endl;

  //Brief: Read Input File With Histogram Binning

  //=======================================
  //     DETECTOR HISTOGRAMS BINNING
  //=======================================
  
  //-------Coincidence-------
  coin_nbins = stod(split(FindString("coin_nbins",    input_HBinFileName.Data())[0], '=')[1]);
  coin_xmin = stod(split(FindString("coin_xmin",      input_HBinFileName.Data())[0], '=')[1]);
  coin_xmax = stod(split(FindString("coin_xmax",      input_HBinFileName.Data())[0], '=')[1]);


  //--HMS DETECTORS--
  hcer_nbins   = stod(split(FindString("hcer_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  hcer_xmin    = stod(split(FindString("hcer_xmin",   input_HBinFileName.Data())[0], '=')[1]);
  hcer_xmax    = stod(split(FindString("hcer_xmax",   input_HBinFileName.Data())[0], '=')[1]);
 
  hcal_nbins   = stod(split(FindString("hcal_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  hcal_xmin    = stod(split(FindString("hcal_xmin",   input_HBinFileName.Data())[0], '=')[1]);
  hcal_xmax    = stod(split(FindString("hcal_xmax",   input_HBinFileName.Data())[0], '=')[1]);

  hbeta_nbins  = stod(split(FindString("hbeta_nbins", input_HBinFileName.Data())[0], '=')[1]);
  hbeta_xmin   = stod(split(FindString("hbeta_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  hbeta_xmax   = stod(split(FindString("hbeta_xmax",  input_HBinFileName.Data())[0], '=')[1]);

  //--SHMS DETECTORS--
  pngcer_nbins = stod(split(FindString("pngcer_nbins",input_HBinFileName.Data())[0], '=')[1]);
  pngcer_xmin  = stod(split(FindString("pngcer_xmin", input_HBinFileName.Data())[0], '=')[1]);
  pngcer_xmax  = stod(split(FindString("pngcer_xmax", input_HBinFileName.Data())[0], '=')[1]);
  
  phgcer_nbins = stod(split(FindString("phgcer_nbins",input_HBinFileName.Data())[0], '=')[1]);
  phgcer_xmin  = stod(split(FindString("phgcer_xmin", input_HBinFileName.Data())[0], '=')[1]);
  phgcer_xmax  = stod(split(FindString("phgcer_xmax", input_HBinFileName.Data())[0], '=')[1]);
 
  pcal_nbins   = stod(split(FindString("pcal_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  pcal_xmin    = stod(split(FindString("pcal_xmin",   input_HBinFileName.Data())[0], '=')[1]);
  pcal_xmax    = stod(split(FindString("pcal_xmax",   input_HBinFileName.Data())[0], '=')[1]);

  pbeta_nbins  = stod(split(FindString("pbeta_nbins", input_HBinFileName.Data())[0], '=')[1]);
  pbeta_xmin   = stod(split(FindString("pbeta_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  pbeta_xmax   = stod(split(FindString("pbeta_xmax",  input_HBinFileName.Data())[0], '=')[1]);

  paero_nbins  = stod(split(FindString("paero_nbins", input_HBinFileName.Data())[0], '=')[1]);
  paero_xmin   = stod(split(FindString("paero_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  paero_xmax   = stod(split(FindString("paero_xmax",  input_HBinFileName.Data())[0], '=')[1]);


  //---------------------------------
  // Kinematics Histograms Binning
  //---------------------------------
  //Primary Kinematics
  the_nbins    	= stod(split(FindString("the_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  the_xmin     	= stod(split(FindString("the_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  the_xmax     	= stod(split(FindString("the_xmax",  input_HBinFileName.Data())[0], '=')[1]);
            				              
  W_nbins      	= stod(split(FindString("W_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  W_xmin       	= stod(split(FindString("W_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  W_xmax       	= stod(split(FindString("W_xmax",  input_HBinFileName.Data())[0], '=')[1]);
               				          
  W2_nbins     	= stod(split(FindString("W2_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  W2_xmin      	= stod(split(FindString("W2_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  W2_xmax      	= stod(split(FindString("W2_xmax",  input_HBinFileName.Data())[0], '=')[1]);
               				          
  Q2_nbins     	= stod(split(FindString("Q2_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  Q2_xmin      	= stod(split(FindString("Q2_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  Q2_xmax      	= stod(split(FindString("Q2_xmax",  input_HBinFileName.Data())[0], '=')[1]);
               				          
  X_nbins      	= stod(split(FindString("X_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  X_xmin       	= stod(split(FindString("X_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  X_xmax       	= stod(split(FindString("X_xmax",  input_HBinFileName.Data())[0], '=')[1]);
  	       				 	  
  nu_nbins     	= stod(split(FindString("nu_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  nu_xmin      	= stod(split(FindString("nu_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  nu_xmax      	= stod(split(FindString("nu_xmax",  input_HBinFileName.Data())[0], '=')[1]);
               				          
  q_nbins      	= stod(split(FindString("q_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  q_xmin       	= stod(split(FindString("q_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  q_xmax       	= stod(split(FindString("q_xmax",  input_HBinFileName.Data())[0], '=')[1]);
               				          
  qx_nbins     	= stod(split(FindString("qx_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  qx_xmin      	= stod(split(FindString("qx_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  qx_xmax      	= stod(split(FindString("qx_xmax",  input_HBinFileName.Data())[0], '=')[1]);
               				          
  qy_nbins     	= stod(split(FindString("qy_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  qy_xmin      	= stod(split(FindString("qy_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  qy_xmax      	= stod(split(FindString("qy_xmax",  input_HBinFileName.Data())[0], '=')[1]);
               				          
  qz_nbins     	= stod(split(FindString("qz_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  qz_xmin      	= stod(split(FindString("qz_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  qz_xmax      	= stod(split(FindString("qz_xmax",  input_HBinFileName.Data())[0], '=')[1]);
               				          
  thq_nbins    	= stod(split(FindString("thq_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  thq_xmin     	= stod(split(FindString("thq_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  thq_xmax     	= stod(split(FindString("thq_xmax",  input_HBinFileName.Data())[0], '=')[1]);
               				          
  phq_nbins    	= stod(split(FindString("phq_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  phq_xmin     	= stod(split(FindString("phq_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  phq_xmax     	= stod(split(FindString("phq_xmax",  input_HBinFileName.Data())[0], '=')[1]);
               				              
  epsilon_nbins	= stod(split(FindString("epsilon_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  epsilon_xmin 	= stod(split(FindString("epsilon_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  epsilon_xmax 	= stod(split(FindString("epsilon_xmax",  input_HBinFileName.Data())[0], '=')[1]);


  //Secondary Kinematics
  Em_nbins     	 = stod(split(FindString("Em_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  Em_xmin      	 = stod(split(FindString("Em_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  Em_xmax      	 = stod(split(FindString("Em_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  	       				  	       
  Pm_nbins     	 = stod(split(FindString("Pm_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  Pm_xmin      	 = stod(split(FindString("Pm_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  Pm_xmax      	 = stod(split(FindString("Pm_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
               				               
  Pmx_lab_nbins	 = stod(split(FindString("Pmx_lab_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  Pmx_lab_xmin 	 = stod(split(FindString("Pmx_lab_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  Pmx_lab_xmax 	 = stod(split(FindString("Pmx_lab_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  	       				  	       
  Pmy_lab_nbins	 = stod(split(FindString("Pmy_lab_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  Pmy_lab_xmin 	 = stod(split(FindString("Pmy_lab_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  Pmy_lab_xmax 	 = stod(split(FindString("Pmy_lab_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  	       				  	       
  Pmz_lab_nbins	 = stod(split(FindString("Pmz_lab_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  Pmz_lab_xmin 	 = stod(split(FindString("Pmz_lab_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  Pmz_lab_xmax 	 = stod(split(FindString("Pmz_lab_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  	       				  	       
  Pmx_q_nbins  	 = stod(split(FindString("Pmx_q_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  Pmx_q_xmin   	 = stod(split(FindString("Pmx_q_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  Pmx_q_xmax   	 = stod(split(FindString("Pmx_q_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  	       				  	       
  Pmy_q_nbins  	 = stod(split(FindString("Pmy_q_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  Pmy_q_xmin   	 = stod(split(FindString("Pmy_q_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  Pmy_q_xmax   	 = stod(split(FindString("Pmy_q_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  	       				  	       
  Pmz_q_nbins  	 = stod(split(FindString("Pmz_q_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  Pmz_q_xmin   	 = stod(split(FindString("Pmz_q_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  Pmz_q_xmax   	 = stod(split(FindString("Pmz_q_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  	       				  	       
  Tx_nbins     	 = stod(split(FindString("Tx_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  Tx_xmin      	 = stod(split(FindString("Tx_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  Tx_xmax      	 = stod(split(FindString("Tx_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  	       				  	       
  Tr_nbins     	 = stod(split(FindString("Tr_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  Tr_xmin      	 = stod(split(FindString("Tr_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  Tr_xmax      	 = stod(split(FindString("Tr_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  	       				  	       
  MM_nbins     	 = stod(split(FindString("MM_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  MM_xmin      	 = stod(split(FindString("MM_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  MM_xmax      	 = stod(split(FindString("MM_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  	       				  	       
  thxq_nbins   	 = stod(split(FindString("thxq_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  thxq_xmin    	 = stod(split(FindString("thxq_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  thxq_xmax    	 = stod(split(FindString("thxq_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  	       				  	       
  thrq_nbins   	 = stod(split(FindString("thrq_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  thrq_xmin    	 = stod(split(FindString("thrq_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  thrq_xmax    	 = stod(split(FindString("thrq_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
               				               
  phxq_nbins   	 = stod(split(FindString("phxq_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  phxq_xmin    	 = stod(split(FindString("phxq_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  phxq_xmax    	 = stod(split(FindString("phxq_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  	       				  	       
  phrq_nbins   	 = stod(split(FindString("phrq_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  phrq_xmin    	 = stod(split(FindString("phrq_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  phrq_xmax    	 = stod(split(FindString("phrq_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  	       				  	       
  Tx_cm_nbins  	 = stod(split(FindString("Tx_cm_nbins", 	input_HBinFileName.Data())[0], '=')[1]);
  Tx_cm_xmin   	 = stod(split(FindString("Tx_cm_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  Tx_cm_xmax   	 = stod(split(FindString("Tx_cm_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  	       				  	       
  Tr_cm_nbins  	 = stod(split(FindString("Tr_cm_nbins", 	input_HBinFileName.Data())[0], '=')[1]);
  Tr_cm_xmin   	 = stod(split(FindString("Tr_cm_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  Tr_cm_xmax   	 = stod(split(FindString("Tr_cm_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  	       				  	       
  thxq_cm_nbins	 = stod(split(FindString("thxq_cm_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  thxq_cm_xmin 	 = stod(split(FindString("thxq_cm_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  thxq_cm_xmax 	 = stod(split(FindString("thxq_cm_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
               				               
  thrq_cm_nbins	 = stod(split(FindString("thrq_cm_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  thrq_cm_xmin 	 = stod(split(FindString("thrq_cm_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  thrq_cm_xmax 	 = stod(split(FindString("thrq_cm_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  	       				  	       
  phxq_cm_nbins	 = stod(split(FindString("phxq_cm_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  phxq_cm_xmin 	 = stod(split(FindString("phxq_cm_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  phxq_cm_xmax 	 = stod(split(FindString("phxq_cm_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  	       				  	       
  phrq_cm_nbins	 = stod(split(FindString("phrq_cm_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  phrq_cm_xmin 	 = stod(split(FindString("phrq_cm_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  phrq_cm_xmax 	 = stod(split(FindString("phrq_cm_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  	       				  	       
  Ttot_cm_nbins	 = stod(split(FindString("Ttot_cm_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  Ttot_cm_xmin 	 = stod(split(FindString("Ttot_cm_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  Ttot_cm_xmax 	 = stod(split(FindString("Ttot_cm_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  	       				  	       
  MandelS_nbins	 = stod(split(FindString("MandelS_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  MandelS_xmin 	 = stod(split(FindString("MandelS_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  MandelS_xmax 	 = stod(split(FindString("MandelS_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  	       				  	       
  MandelT_nbins	 = stod(split(FindString("MandelT_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  MandelT_xmin 	 = stod(split(FindString("MandelT_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  MandelT_xmax 	 = stod(split(FindString("MandelT_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  	       				  	       
  MandelU_nbins	 = stod(split(FindString("MandelU_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  MandelU_xmin 	 = stod(split(FindString("MandelU_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  MandelU_xmax 	 = stod(split(FindString("MandelU_xmax",  	input_HBinFileName.Data())[0], '=')[1]);


    
  //Kinematics Defined in HCANA (which are not in primary/secondary modules)
  kf_nbins	= stod(split(FindString("kf_nbins",  	input_HBinFileName.Data())[0], '=')[1]);   //final electron momentum
  kf_xmin	= stod(split(FindString("kf_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  kf_xmax	= stod(split(FindString("kf_xmax",  	input_HBinFileName.Data())[0], '=')[1]);

  Pf_nbins	= stod(split(FindString("Pf_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  Pf_xmin	= stod(split(FindString("Pf_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  Pf_xmax	= stod(split(FindString("Pf_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  
  //Additional Kinematics
  thx_nbins	= stod(split(FindString("thx_nbins",  	input_HBinFileName.Data())[0], '=')[1]);  //proton(hadron) angle
  thx_xmin 	= stod(split(FindString("thx_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  thx_xmax 	= stod(split(FindString("thx_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
           				          
  MM2_nbins	= stod(split(FindString("MM2_nbins",  	input_HBinFileName.Data())[0], '=')[1]);  
  MM2_xmin 	= stod(split(FindString("MM2_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  MM2_xmax 	= stod(split(FindString("MM2_xmax",  	input_HBinFileName.Data())[0], '=')[1]);

  //---------------------------------
  // Acceptance Histograms Binning
  //---------------------------------

  //----Electron Arm Focal Plane-----
  exfp_nbins 	= stod(split(FindString("exfp_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  exfp_xmin  	= stod(split(FindString("exfp_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  exfp_xmax  	= stod(split(FindString("exfp_xmax",  input_HBinFileName.Data())[0], '=')[1]);
             				            
  expfp_nbins	= stod(split(FindString("expfp_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  expfp_xmin 	= stod(split(FindString("expfp_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  expfp_xmax 	= stod(split(FindString("expfp_xmax",  input_HBinFileName.Data())[0], '=')[1]);
  	     				 	    
  eyfp_nbins 	= stod(split(FindString("eyfp_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  eyfp_xmin  	= stod(split(FindString("eyfp_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  eyfp_xmax  	= stod(split(FindString("eyfp_xmax",  input_HBinFileName.Data())[0], '=')[1]);
             				            
  eypfp_nbins	= stod(split(FindString("eypfp_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  eypfp_xmin 	= stod(split(FindString("eypfp_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  eypfp_xmax 	= stod(split(FindString("eypfp_xmax",  input_HBinFileName.Data())[0], '=')[1]);

  //----Electron Arm Reconstructed-----
  eytar_nbins 	= stod(split(FindString("eytar_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  eytar_xmin  	= stod(split(FindString("eytar_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  eytar_xmax  	= stod(split(FindString("eytar_xmax",  input_HBinFileName.Data())[0], '=')[1]);
              				             
  eyptar_nbins	= stod(split(FindString("eyptar_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  eyptar_xmin 	= stod(split(FindString("eyptar_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  eyptar_xmax 	= stod(split(FindString("eyptar_xmax",  input_HBinFileName.Data())[0], '=')[1]);
  	      				 	     
  exptar_nbins	= stod(split(FindString("exptar_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  exptar_xmin 	= stod(split(FindString("exptar_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  exptar_xmax 	= stod(split(FindString("exptar_xmax",  input_HBinFileName.Data())[0], '=')[1]);
  	      				 	     
  edelta_nbins	= stod(split(FindString("edelta_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  edelta_xmin 	= stod(split(FindString("edelta_xmin",  input_HBinFileName.Data())[0], '=')[1]);  
  edelta_xmax 	= stod(split(FindString("edelta_xmax",  input_HBinFileName.Data())[0], '=')[1]);   

  
  //----Hadron Arm Focal Plane-----
  hxfp_nbins 	= stod(split(FindString("hxfp_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  hxfp_xmin  	= stod(split(FindString("hxfp_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  hxfp_xmax  	= stod(split(FindString("hxfp_xmax",  input_HBinFileName.Data())[0], '=')[1]);
             				            
  hxpfp_nbins	= stod(split(FindString("hxpfp_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  hxpfp_xmin 	= stod(split(FindString("hxpfp_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  hxpfp_xmax 	= stod(split(FindString("hxpfp_xmax",  input_HBinFileName.Data())[0], '=')[1]);
  	     				 	    
  hyfp_nbins 	= stod(split(FindString("hyfp_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  hyfp_xmin  	= stod(split(FindString("hyfp_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  hyfp_xmax  	= stod(split(FindString("hyfp_xmax",  input_HBinFileName.Data())[0], '=')[1]);
             				            
  hypfp_nbins	= stod(split(FindString("hypfp_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  hypfp_xmin 	= stod(split(FindString("hypfp_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  hypfp_xmax 	= stod(split(FindString("hypfp_xmax",  input_HBinFileName.Data())[0], '=')[1]);

  //----Hadron Arm Reconstructed-----
  hytar_nbins 	= stod(split(FindString("hytar_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  hytar_xmin  	= stod(split(FindString("hytar_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  hytar_xmax  	= stod(split(FindString("hytar_xmax",  input_HBinFileName.Data())[0], '=')[1]);
              				             
  hyptar_nbins	= stod(split(FindString("hyptar_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  hyptar_xmin 	= stod(split(FindString("hyptar_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  hyptar_xmax 	= stod(split(FindString("hyptar_xmax",  input_HBinFileName.Data())[0], '=')[1]);
  	      				 	     
  hxptar_nbins	= stod(split(FindString("hxptar_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  hxptar_xmin 	= stod(split(FindString("hxptar_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  hxptar_xmax 	= stod(split(FindString("hxptar_xmax",  input_HBinFileName.Data())[0], '=')[1]);
  	      				 	     
  hdelta_nbins	= stod(split(FindString("hdelta_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  hdelta_xmin 	= stod(split(FindString("hdelta_xmin",  input_HBinFileName.Data())[0], '=')[1]);  
  hdelta_xmax 	= stod(split(FindString("hdelta_xmax",  input_HBinFileName.Data())[0], '=')[1]);   


  //----Target Quantities----
  //(Use same binning for hadron/electron reconstructed at target)
  tarx_nbins	= stod(split(FindString("tarx_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  tarx_xmin 	= stod(split(FindString("tarx_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  tarx_xmax 	= stod(split(FindString("tarx_xmax",  input_HBinFileName.Data())[0], '=')[1]);
            				           
  tary_nbins	= stod(split(FindString("tary_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  tary_xmin 	= stod(split(FindString("tary_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  tary_xmax 	= stod(split(FindString("tary_xmax",  input_HBinFileName.Data())[0], '=')[1]);

  tarz_nbins	= stod(split(FindString("tarz_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  tarz_xmin 	= stod(split(FindString("tarz_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  tarz_xmax 	= stod(split(FindString("tarz_xmax",  input_HBinFileName.Data())[0], '=')[1]);

  ztar_diff_nbins	= stod(split(FindString("ztar_diff_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  ztar_diff_xmin 	= stod(split(FindString("ztar_diff_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  ztar_diff_xmax 	= stod(split(FindString("ztar_diff_xmax",  input_HBinFileName.Data())[0], '=')[1]);

  
  //----Collimator Quantities----
  hXColl_nbins	= stod(split(FindString("hXColl_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  hXColl_xmin 	= stod(split(FindString("hXColl_xmin",  input_HBinFileName.Data())[0], '=')[1]);  
  hXColl_xmax 	= stod(split(FindString("hXColl_xmax",  input_HBinFileName.Data())[0], '=')[1]);   
  	      				 	     
  hYColl_nbins	= stod(split(FindString("hYColl_nbins",  input_HBinFileName.Data())[0], '=')[1]);                                           
  hYColl_xmin 	= stod(split(FindString("hYColl_xmin",  input_HBinFileName.Data())[0], '=')[1]);                                                     
  hYColl_xmax 	= stod(split(FindString("hYColl_xmax",  input_HBinFileName.Data())[0], '=')[1]);
  	      				 	     
  eXColl_nbins	= stod(split(FindString("eXColl_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  eXColl_xmin 	= stod(split(FindString("eXColl_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  eXColl_xmax 	= stod(split(FindString("eXColl_xmax",  input_HBinFileName.Data())[0], '=')[1]);
  	      				 	     
  eYColl_nbins	= stod(split(FindString("eYColl_nbins",  input_HBinFileName.Data())[0], '=')[1]);      
  eYColl_xmin 	= stod(split(FindString("eYColl_xmin",  input_HBinFileName.Data())[0], '=')[1]);                                                                      
  eYColl_xmax 	= stod(split(FindString("eYColl_xmax",  input_HBinFileName.Data())[0], '=')[1]);

  

}

//_______________________________________________________________________________
void baseAnalyzer::CreateHist()
{

  cout << "Calling Base CreateHist()  " << endl;
  //Method to create Histograms

  //Create TLists to store categorical histograms
  pid_HList  = new TList();
  kin_HList  = new TList();
  accp_HList = new TList();

  
  // Dummy histograms to store the ith and cumulative histograms (See CombineHistos() Method)
  //1D
  h_total = new TH1F();         //dummy histo to store hsitogram sum
  h_i     = new TH1F();         //dummy histo to store ith histogram from list
  //2D
  h2_total = new TH2F();        //dummy histo to store hsitogram sum
  h2_i     = new TH2F();       //dummy histo to store ith histogram from list

  
  //--------------------------------------------------------------------
  //---------HISTOGRAM CATEGORY: Particle Identification (PID)----------
  //--------------------------------------------------------------------
  
  //Coincidence HISTOS
  H_ep_ctime   = new TH1F("H_ep_ctime", "ep Coincidence Time; ep Coincidence Time [ns]; Counts / mC", coin_nbins, coin_xmin, coin_xmax);
  H_ep_ctime->Sumw2(); //Apply sum of weight squared to this histogram ABOVE.
  H_ep_ctime->SetDefaultSumw2(kTRUE);  //Generalize sum weights squared to all histograms  (ROOT 6 has this by default. ROOT 5 does NOT)
  
  H_eK_ctime   = new TH1F("H_eK_ctime", "eK Coincidence Time; eK Coincidence Time [ns]; Counts / mC", coin_nbins, coin_xmin, coin_xmax);
  H_ePi_ctime  = new TH1F("H_ePi_ctime", "e#pi Coincidence Time; e#pi Coincidence Time [ns]; Counts / mC", coin_nbins, coin_xmin, coin_xmax);
    
  //HMS DETECTORS HISTOS
  H_hCerNpeSum      = new TH1F("H_hCerNpeSum", "HMS Cherenkov NPE Sum; Cherenkov NPE Sum; Counts / mC", hcer_nbins, hcer_xmin, hcer_xmax);
  H_hCalEtotNorm    = new TH1F("H_hCalEtotNorm", "HMS Calorimeter Normalized Total Energy; E_{tot} / P_{cent}; Counts / mC", hcal_nbins, hcal_xmin, hcal_xmax);
  H_hCalEtotTrkNorm = new TH1F("H_hCalEtotTrkNorm", "HMS Calorimeter Total Normalized Track Energy; E_{tot} / P_{trk}; Counts / mC", hcal_nbins, hcal_xmin, hcal_xmax);
  H_hHodBetaNtrk    = new TH1F("H_hHodBetaNtrk", "HMS Hodo #beta (no track); #beta (no track); Counts / mC", hbeta_nbins, hbeta_xmin, hbeta_xmax);
  H_hHodBetaTrk     = new TH1F("H_hHodBetaTrk", "HMS Hodo #beta (golden track); #beta (golden track); Counts / mC", hbeta_nbins, hbeta_xmin, hbeta_xmax);
  
  //SHMS DETECTORS HISTOS
  H_pNGCerNpeSum    = new TH1F("H_pNGCerNpeSum", "SHMS Noble Gas Cherenkov NPE Sum; Cherenkov NPE Sum; Counts / mC ", pngcer_nbins, pngcer_xmin, pngcer_xmax);
  H_pHGCerNpeSum    = new TH1F("H_pHGCerNpeSum", "SHMS Heavy Gas Cherenkov NPE Sum; Cherenkov NPE Sum; Counts / mC ", phgcer_nbins, phgcer_xmin, phgcer_xmax);
  H_pAeroNpeSum     = new TH1F("H_pAeroNpeSum", "SHMS Aerogel NPE Sum; Aerogel NPE Sum; Counts / mC ", paero_nbins, paero_xmin, paero_xmax);
  H_pCalEtotNorm    = new TH1F("H_pCalEtotNorm", "SHMS Calorimeter Normalized Total Energy; E_{tot} / P_{cent}; Counts / mC", pcal_nbins, pcal_xmin, pcal_xmax);
  H_pCalEtotTrkNorm = new TH1F("H_pCalEtotTrkNorm", "SHMS Calorimeter Total Normalized Track Energy; E_{tot} / P_{trk}; Counts / mC", pcal_nbins, pcal_xmin, pcal_xmax);
  H_pHodBetaNtrk    = new TH1F("H_pBetaNtrk", "SHMS Hodo #beta (no track); #beta (no track); Counts / mC", pbeta_nbins, pbeta_xmin, pbeta_xmax);
  H_pHodBetaTrk     = new TH1F("H_pBetaTrk", "SHMS Hodo #beta (golden track); #beta (golden track); Counts / mC", pbeta_nbins, pbeta_xmin, pbeta_xmax);

  //HMS 2D PID               
  H_hcal_vs_hcer     = new TH2F("H_hcal_vs_hcer", "HMS: Calorimeter vs. Cherenkov; Calorimeter E_{tot}/P_{trk}; Cherenkov NPE Sum", hcal_nbins, hcal_xmin, hcal_xmax, hcer_nbins, hcer_xmin, hcer_xmax);     
  		     	     
  //SHMS 2D PID	     	     
  H_pcal_vs_phgcer    = new TH2F("H_pcal_vs_phgcer", "SHMS: Heavy Gas Cherenkov (HGC) vs. Calorimeter; Calorimeter E_{tot}/P_{trk}; HGC NPE Sum", pcal_nbins, pcal_xmin, pcal_xmax, phgcer_nbins, phgcer_xmin, phgcer_xmax);        
  H_pcal_vs_pngcer    = new TH2F("H_pcal_vs_pngcer", "SHMS: Noble Gas Cherenkov (NGC) vs. Calorimeter; Calorimeter E_{tot}/P_{trk}; NGC NPE Sum", pcal_nbins, pcal_xmin, pcal_xmax, pngcer_nbins, pngcer_xmin, pngcer_xmax);   
  H_pcal_vs_paero     = new TH2F("H_pcal_vs_paero",  "SHMS: Aerogel Cherenkov (AER) vs. Calorimeter;   Calorimeter E_{tot}/P_{trk}; AER NPE Sum", pcal_nbins, pcal_xmin, pcal_xmax, paero_nbins,  paero_xmin,  paero_xmax);    
  H_paero_vs_phgcer   = new TH2F("H_paero_vs_phgcer","SHMS: Heavy Gas Cherenkov (HGC) vs. Aerogel (AER); AER NPE Sum; HGC NPE Sum", paero_nbins,  paero_xmin,  paero_xmax, phgcer_nbins, phgcer_xmin, phgcer_xmax);      
  H_paero_vs_pngcer   = new TH2F("H_paero_vs_pngcer","SHMS: Noble Gas Cherenkov (NGC) vs. Aerogel (AER); AER NPE Sum; NGC NPE Sum", paero_nbins,  paero_xmin,  paero_xmax, pngcer_nbins, pngcer_xmin, pngcer_xmax);   
  H_pngcer_vs_phgcer  = new TH2F("H_pngcer_vs_phgcer","SHMS: Heavy Gas Cherenkov (HGC) vs. Noble Gas Cherenkov (NGC); NGC NPE Sum; HGC NPE Sum", pngcer_nbins, pngcer_xmin, pngcer_xmax, phgcer_nbins, phgcer_xmin, phgcer_xmax); 
  
  
  //Add PID Histos to TList
  pid_HList->Add(H_ep_ctime);
  pid_HList->Add(H_eK_ctime);
  pid_HList->Add(H_ePi_ctime);
  pid_HList->Add(H_hCerNpeSum);
  pid_HList->Add(H_hCalEtotNorm);
  pid_HList->Add(H_hCalEtotTrkNorm);
  pid_HList->Add(H_hHodBetaNtrk);
  pid_HList->Add(H_hHodBetaTrk);
  pid_HList->Add(H_pNGCerNpeSum);
  pid_HList->Add(H_pHGCerNpeSum);
  pid_HList->Add(H_pAeroNpeSum);
  pid_HList->Add(H_pCalEtotNorm);
  pid_HList->Add(H_pCalEtotTrkNorm);
  pid_HList->Add(H_pHodBetaNtrk);
  pid_HList->Add(H_pHodBetaTrk);
  //Add 2D PID
  pid_HList->Add(H_hcal_vs_hcer);  
  pid_HList->Add(H_pcal_vs_phgcer);
  pid_HList->Add(H_pcal_vs_pngcer);
  pid_HList->Add(H_pcal_vs_paero);
  pid_HList->Add(H_paero_vs_phgcer);
  pid_HList->Add(H_paero_vs_pngcer);
  pid_HList->Add(H_pngcer_vs_phgcer);

  
  
  //--------------------------------------------------------
  //---------HISTOGRAM CATEGORY: Kinematics  (KIN)----------
  //--------------------------------------------------------

  //Primary (electron) Kinematics (14 histos)
  H_the     = new TH1F("H_the", "Electron Scattering Angle, #theta_{e}", the_nbins, the_xmin, the_xmax);
  H_kf      = new TH1F("H_kf", "Final e^{-} Momentum", kf_nbins, kf_xmin, kf_xmax);
  H_W       = new TH1F("H_W", "Invariant Mass, W", W_nbins, W_xmin, W_xmax); 
  H_W2      = new TH1F("H_W2", "Invariant Mass, W^{2}", W2_nbins, W2_xmin, W2_xmax); 
  H_Q2      = new TH1F("H_Q2","4-Momentum Transfer, Q^{2}", Q2_nbins, Q2_xmin, Q2_xmax); 
  H_xbj     = new TH1F("H_xbj", "x-Bjorken", X_nbins, X_xmin, X_xmax);  
  H_nu      = new TH1F("H_nu","Energy Transfer, #nu", nu_nbins, nu_xmin, nu_xmax); 
  H_q       = new TH1F("H_q", "3-Momentum Transfer, |#vec{q}|", q_nbins, q_xmin, q_xmax);
  H_qx      = new TH1F("H_qx", "|#vec{q}_{x}|", qx_nbins, qx_xmin, qx_xmax);
  H_qy      = new TH1F("H_qy", "|#vec{q}_{y}|", qy_nbins, qy_xmin, qy_xmax);
  H_qz      = new TH1F("H_qz", "|#vec{q}_{z}|", qz_nbins, qz_xmin, qz_xmax);
  H_thq     = new TH1F("H_thq", "In-Plane Angle w.r.t +z(lab), #theta_{q}", thq_nbins, thq_xmin, thq_xmax); 
  H_phq     = new TH1F("H_phq", "Out-of-Plane Angle w.r.t +z(lab), #phi_{q}", phq_nbins, phq_xmin, phq_xmax); 
  H_epsilon = new TH1F("H_epsilon", "Virtual Photon (#gamma) Polarization Factor" , epsilon_nbins, epsilon_xmin, epsilon_xmax); 

  //Secondary (Hadron) Kinematics (recoil and missing are used interchageably) ()
  H_Em      = new TH1F("H_Emiss","Missing Energy", Em_nbins, Em_xmin, Em_xmax);   
  H_Em_nuc  = new TH1F("H_Em_nuc","Nuclear Missing Energy", Em_nbins, Em_xmin, Em_xmax); 
  H_Pm      = new TH1F("H_Pm","Missing Momentum, P_{miss}", Pm_nbins, Pm_xmin, Pm_xmax); 
  H_Pmx_lab = new TH1F("H_Pmx_Lab","P_{miss, x} (Lab)", Pmx_lab_nbins, Pmx_lab_xmin, Pmx_lab_xmax);         
  H_Pmy_lab = new TH1F("H_Pmy_Lab","P_{miss, y} (Lab)", Pmy_lab_nbins, Pmy_lab_xmin, Pmy_lab_xmax);    
  H_Pmz_lab = new TH1F("H_Pmz_Lab","P_{miss, z} (Lab)", Pmz_lab_nbins, Pmz_lab_xmin, Pmz_lab_xmax);  
  H_Pmx_q   = new TH1F("H_Pmx_q","P_{miss, xq} (w.r.t #vec{q}) ", Pmx_q_nbins, Pmx_q_xmin, Pmx_q_xmax);   
  H_Pmy_q   = new TH1F("H_Pmy_q","P_{miss, yq} (w.r.t #vec{q}) ", Pmy_q_nbins, Pmy_q_xmin, Pmy_q_xmax); 
  H_Pmz_q   = new TH1F("H_Pmz_q","P_{miss, zq} (along #vec{q}) ", Pmz_q_nbins, Pmz_q_xmin, Pmz_q_xmax); 
  H_Tx      = new TH1F("H_Tx", "Kinetic Energy, T_{x} (detected)", Tx_nbins, Tx_xmin, Tx_xmax);     
  H_Tr      = new TH1F("H_Tr", "Kinetic Energy, T_{r} (recoil)",   Tr_nbins, Tr_xmin, Tr_xmax);  
  H_MM      = new TH1F("H_MM","Missing Mass, M_{miss}", MM_nbins, MM_xmin, MM_xmax);        
  H_MM2     = new TH1F("H_MM2","Missing Mass Squared, M^{2}_{miss}", MM2_nbins, MM2_xmin, MM2_xmax); 
  H_thx     = new TH1F("H_thx", "Hadron Scattering Angle (detected), #theta_{x}", thx_nbins, thx_xmin, thx_xmax);
  H_Pf      = new TH1F("H_Pf", "Final Hadron Momentum (detected), p_{f}", Pf_nbins, Pf_xmin, Pf_xmax);
  H_thxq    = new TH1F("H_thxq", "In-Plane Angle, #theta_{xq}", thxq_nbins, thxq_xmin, thxq_xmax);
  H_thrq    = new TH1F("H_thrq", "In-Plane Angle, #theta_{rq}", thrq_nbins, thrq_xmin, thrq_xmax);
  H_phxq    = new TH1F("H_phxq", "Out-of-Plane Angle, #phi_{xq}", phxq_nbins, phxq_xmin, phxq_xmax);
  H_phrq    = new TH1F("H_phrq", "Out-of-Plane Angle, #phi_{rq}", phrq_nbins, phrq_xmin, phrq_xmax);
  H_Tx_cm   = new TH1F("H_Tx_cm", "Kinetic Energy, T_{x, cm} (detected)", Tx_cm_nbins, Tx_cm_xmin, Tx_cm_xmax);     
  H_Tr_cm   = new TH1F("H_Tr_cm", "Kinetic Energy, T_{r, cm} (recoil)",   Tr_cm_nbins, Tr_cm_xmin, Tr_cm_xmax);  
  H_thxq_cm = new TH1F("H_thxq_cm", "In-Plane Angle, #theta_{xq, cm}", thxq_cm_nbins, thxq_cm_xmin, thxq_cm_xmax);
  H_thrq_cm = new TH1F("H_thrq_cm", "In-Plane Angle, #theta_{rq, cm}", thrq_cm_nbins, thrq_cm_xmin, thrq_cm_xmax);
  H_phxq_cm = new TH1F("H_phxq_cm", "Out-of-Plane Angle, #phi_{xq, cm}", phxq_cm_nbins, phxq_cm_xmin, phxq_cm_xmax);
  H_phrq_cm = new TH1F("H_phrq_cm", "Out-of-Plane Angle, #phi_{rq, cm}", phrq_cm_nbins, phrq_cm_xmin, phrq_cm_xmax);
  H_Ttot_cm = new TH1F("H_Ttot_cm", "Total CM Kinetic Energy, T_{tot,cm}", Ttot_cm_nbins, Ttot_cm_xmin, Ttot_cm_xmax);
  H_MandelS = new TH1F("H_MandelS", "s-Mandelstam", MandelS_nbins, MandelS_xmin, MandelS_xmax);     
  H_MandelT = new TH1F("H_MandelT", "t-Mandelstam", MandelT_nbins, MandelT_xmin, MandelT_xmax);
  H_MandelU = new TH1F("H_MandelU", "u-Mandelstam", MandelU_nbins, MandelU_xmin, MandelU_xmax);     
  
  // (Cosine, Sine) Histos of detected AND recoil angles (range is fixed at: -1, 1)
  //LAB FRAME
  H_cth_xq = new TH1F("H_cth_xq", "cos(#theta_{xq})", thxq_nbins, -1, 1);
  H_cth_rq = new TH1F("H_cth_rq", "cos(#theta_{rq})", thrq_nbins, -1, 1);
  H_sth_xq = new TH1F("H_sth_xq", "sin(#theta_{xq})", thxq_nbins, -1, 1);
  H_sth_rq = new TH1F("H_sth_rq", "sin(#theta_{rq})", thrq_nbins, -1, 1);
  H_cphi_xq = new TH1F("H_cphi_xq", "cos(#phi_{xq})", phxq_nbins, -1, 1);
  H_cphi_rq = new TH1F("H_cphi_rq", "cos(#phi_{rq})", phrq_nbins, -1, 1);
  H_sphi_xq = new TH1F("H_sphi_xq", "sin(#phi_{xq})", phxq_nbins, -1, 1);
  H_sphi_rq = new TH1F("H_sphi_rq", "sin(#phi_{rq})", phrq_nbins, -1, 1);
  //CM FRAME
  H_cth_xq_cm = new TH1F("H_cth_xq_cm", "cos(#theta_{xq,cm})", thxq_cm_nbins, -1, 1);
  H_cth_rq_cm = new TH1F("H_cth_rq_cm", "cos(#theta_{rq,cm})", thrq_cm_nbins, -1, 1);
  H_sth_xq_cm = new TH1F("H_sth_xq_cm", "sin(#theta_{xq,cm})", thxq_cm_nbins, -1, 1);
  H_sth_rq_cm = new TH1F("H_sth_rq_cm", "sin(#theta_{rq,cm})", thrq_cm_nbins, -1, 1);
  H_cphi_xq_cm = new TH1F("H_cphi_xq_cm", "cos(#phi_{xq,cm})", phxq_cm_nbins, -1, 1);
  H_cphi_rq_cm = new TH1F("H_cphi_rq_cm", "cos(#phi_{rq,cm})", phrq_cm_nbins, -1, 1);
  H_sphi_xq_cm = new TH1F("H_sphi_xq_cm", "sin(#phi_{xq,cm})", phxq_cm_nbins, -1, 1);
  H_sphi_rq_cm = new TH1F("H_sphi_rq_cm", "sin(#phi_{rq,cm})", phrq_cm_nbins, -1, 1);
  
  //Add Kin Histos to TList

  //Add Primary Kin Histos
  kin_HList->Add( H_the    );
  kin_HList->Add( H_kf     );
  kin_HList->Add( H_W      );
  kin_HList->Add( H_W2     );
  kin_HList->Add( H_Q2     );
  kin_HList->Add( H_xbj    );
  kin_HList->Add( H_nu     );
  kin_HList->Add( H_q      );
  kin_HList->Add( H_qx     );
  kin_HList->Add( H_qy     );
  kin_HList->Add( H_qz     );
  kin_HList->Add( H_thq    );
  kin_HList->Add( H_phq    );
  kin_HList->Add( H_epsilon); 

  //Add Secondary Kin Histos
  kin_HList->Add( H_Em       );
  kin_HList->Add( H_Em_nuc   );
  kin_HList->Add( H_Pm       );
  kin_HList->Add( H_Pmx_lab  );
  kin_HList->Add( H_Pmy_lab  );
  kin_HList->Add( H_Pmz_lab  );
  kin_HList->Add( H_Pmx_q    );
  kin_HList->Add( H_Pmy_q    );
  kin_HList->Add( H_Pmz_q    );
  kin_HList->Add( H_Tx       );
  kin_HList->Add( H_Tr       );
  kin_HList->Add( H_MM       );
  kin_HList->Add( H_MM2      );
  kin_HList->Add( H_thx      );
  kin_HList->Add( H_Pf       );
  kin_HList->Add( H_thxq     );
  kin_HList->Add( H_thrq     );
  kin_HList->Add( H_phxq     );
  kin_HList->Add( H_phrq     );
  kin_HList->Add( H_Tx_cm    );
  kin_HList->Add( H_Tr_cm    );
  kin_HList->Add( H_thxq_cm  );
  kin_HList->Add( H_thrq_cm  );
  kin_HList->Add( H_phxq_cm  );
  kin_HList->Add( H_phrq_cm  );
  kin_HList->Add( H_Ttot_cm  );
  kin_HList->Add( H_MandelS  );
  kin_HList->Add( H_MandelT  );
  kin_HList->Add( H_MandelU  );

  //Add (cosine, sine) of angles relative to q
  //LAB FRAME
  kin_HList->Add( H_cth_xq  );
  kin_HList->Add( H_cth_rq  );
  kin_HList->Add( H_sth_xq  );
  kin_HList->Add( H_sth_rq  );
  kin_HList->Add( H_cphi_xq );
  kin_HList->Add( H_cphi_rq );
  kin_HList->Add( H_sphi_xq );
  kin_HList->Add( H_sphi_rq );
  //CM FRAME
  kin_HList->Add( H_cth_xq_cm  );
  kin_HList->Add( H_cth_rq_cm  );
  kin_HList->Add( H_sth_xq_cm  );
  kin_HList->Add( H_sth_rq_cm  );
  kin_HList->Add( H_cphi_xq_cm );
  kin_HList->Add( H_cphi_rq_cm );
  kin_HList->Add( H_sphi_xq_cm );
  kin_HList->Add( H_sphi_rq_cm );


  //----------------------------------------------------------------------
  //---------HISTOGRAM CATEGORY: Spectrometer Acceptance  (ACCP)----------
  //----------------------------------------------------------------------


  //Electron Arm Focal Plane Quantities
  H_exfp = new TH1F("H_exfp", Form("%s X_{fp}; X_{fp} [cm]; Counts / mC", e_arm_name.Data()), exfp_nbins, exfp_xmin, exfp_xmax);
  H_eyfp = new TH1F("H_eyfp", Form("%s Y_{fp}; Y_{fp} [cm]; Counts / mC", e_arm_name.Data()), eyfp_nbins, eyfp_xmin, eyfp_xmax);
  H_expfp = new TH1F("H_expfp", Form("%s X'_{fp}; X'_{fp} [rad]; Counts / mC", e_arm_name.Data()), expfp_nbins, expfp_xmin, expfp_xmax);
  H_eypfp = new TH1F("H_eypfp", Form("%s Y'_{fp}; Y'_{fp} [rad]; Counts / mC", e_arm_name.Data()), eypfp_nbins, eypfp_xmin, eypfp_xmax);
  
  //Electron Arm Reconstructed Quantities 
  H_eytar = new TH1F("H_eytar", Form("%s Y_{tar}; Y_{tar} [cm]; Counts / mC", e_arm_name.Data()), eytar_nbins, eytar_xmin, eytar_xmax);
  H_exptar = new TH1F("H_exptar", Form("%s X'_{tar}; X'_{tar} [rad]; Counts / mC", e_arm_name.Data()), exptar_nbins, exptar_xmin, exptar_xmax);
  H_eyptar = new TH1F("H_eyptar", Form("%s Y'_{tar}; Y'_{tar} [rad]; Counts / mC", e_arm_name.Data()), eyptar_nbins, eyptar_xmin, eyptar_xmax);
  H_edelta = new TH1F("H_edelta", Form("%s Momentum Acceptance, #delta; #delta [%%]; Counts / mC", e_arm_name.Data()), edelta_nbins, edelta_xmin, edelta_xmax);
  
  //Hadron arm Focal Plane Quantities
  H_hxfp = new TH1F("H_hxfp", Form("%s  X_{fp}; X_{fp} [cm]; Counts / mC", h_arm_name.Data()), hxfp_nbins, hxfp_xmin, hxfp_xmax);
  H_hyfp = new TH1F("H_hyfp", Form("%s  Y_{fp}; Y_{fp} [cm]; Counts / mC", h_arm_name.Data()), hyfp_nbins, hyfp_xmin, hyfp_xmax);
  H_hxpfp = new TH1F("H_hxpfp", Form("%s  X'_{fp}; X'_{fp} [rad]; Counts / mC", h_arm_name.Data()), hxpfp_nbins, hxpfp_xmin, hxpfp_xmax );
  H_hypfp = new TH1F("H_hypfp", Form("%s  Y'_{fp}; Y'_{fp} [rad]; Counts / mC", h_arm_name.Data()), hypfp_nbins, hypfp_xmin, hypfp_xmax);

  //Hadron arm Reconstructed Quantities 
  H_hytar = new TH1F("H_hytar", Form("%s  Y_{tar}; Y_{tar} [cm]; Counts / mC", h_arm_name.Data()), hytar_nbins, hytar_xmin, hytar_xmax);
  H_hxptar = new TH1F("H_hxptar", Form("%s  X'_{tar}; X'_{tar} [rad]; Counts / mC", h_arm_name.Data()), hxptar_nbins, hxptar_xmin, hxptar_xmax);
  H_hyptar = new TH1F("H_hyptar", Form("%s  Y'_{tar}; Y'_{tar} [rad]; Counts / mC", h_arm_name.Data()), hyptar_nbins, hyptar_xmin, hyptar_xmax );
  H_hdelta = new TH1F("H_hdelta", Form("%s  Momentum Acceptance, #delta; #delta [%%]; Counts / mC", h_arm_name.Data()), hdelta_nbins, hdelta_xmin, hdelta_xmax);
  

  //Target Reconstruction (Hall Coord. System) 
  H_htar_x = new TH1F("H_htar_x", Form("%s x-Target (Lab); x-Target [cm]; Counts / mC", h_arm_name.Data()), tarx_nbins, tarx_xmin, tarx_xmax);
  H_htar_y = new TH1F("H_htar_y", Form("%s y_Target (Lab); y-Target [cm]; Counts / mC", h_arm_name.Data()), tary_nbins, tary_xmin, tary_xmax);
  H_htar_z = new TH1F("H_htar_z", Form("%s z_Target (Lab); z-Target [cm]; Counts / mC", h_arm_name.Data()), tarz_nbins, tarz_xmin, tarz_xmax);
  H_etar_x = new TH1F("H_etar_x", Form("%s x-Target (Lab); x-Target [cm]; Counts / mC", e_arm_name.Data()), tarx_nbins, tarx_xmin, tarx_xmax);
  H_etar_y = new TH1F("H_etar_y", Form("%s y-Target (Lab); y-Target [cm]; Counts / mC", e_arm_name.Data()), tary_nbins, tary_xmin, tary_xmax);
  H_etar_z = new TH1F("H_etar_z", Form("%s z-Target (Lab); z-Target [cm]; Counts / mC", e_arm_name.Data()), tarz_nbins, tarz_xmin, tarz_xmax);

  //difference in reaction vertex z (user-defined)
  H_ztar_diff = new TH1F("H_ztar_diff", "Ztar Difference; z-Target Difference [cm]; Counts / mC", ztar_diff_nbins, ztar_diff_xmin, ztar_diff_xmax);

  //HMS / SHMS Collimator
  H_hXColl = new TH1F("H_hXColl", Form("%s X Collimator; X-Collimator [cm]; Counts / mC", h_arm_name.Data()), hXColl_nbins, hXColl_xmin, hXColl_xmax);
  H_hYColl = new TH1F("H_hYColl", Form("%s Y Collimator; Y-Collimator [cm]; Counts / mC", h_arm_name.Data()), hYColl_nbins, hYColl_xmin, hYColl_xmax); 
  H_eXColl = new TH1F("H_eXColl", Form("%s X Collimator; X-Collimator [cm]; Counts / mC", e_arm_name.Data()), eXColl_nbins, eXColl_xmin, eXColl_xmax);                                                                             
  H_eYColl = new TH1F("H_eYColl", Form("%s Y Collimator; Y-Collimator [cm]; Counts / mC", e_arm_name.Data()), eYColl_nbins, eYColl_xmin, eYColl_xmax);        

  //2D Collimator Histos
  H_hXColl_vs_hYColl = new TH2F("H_hXColl_vs_hYColl", Form("%s Collimator; %s Y-Collimator [cm]; %s X-Collimator [cm]", h_arm_name.Data(), h_arm_name.Data(), h_arm_name.Data()), hYColl_nbins, hYColl_xmin, hYColl_xmax,  hXColl_nbins, hXColl_xmin, hXColl_xmax);
  H_eXColl_vs_eYColl = new TH2F("H_eXColl_vs_eYColl", Form("%s Collimator; %s Y-Collimator [cm]; %s X-Collimator [cm]", e_arm_name.Data(), e_arm_name.Data(), e_arm_name.Data()), eYColl_nbins, eYColl_xmin, eYColl_xmax, eXColl_nbins, eXColl_xmin, eXColl_xmax); 
  
  //2D Hour Glass Histos
  H_hxfp_vs_hyfp  = new TH2F("H_hxfp_vs_hyfp", Form("%s  X_{fp} vs. Y_{fp}; Y_{fp} [cm]; X_{fp} [cm]", h_arm_name.Data()),  hyfp_nbins, hyfp_xmin, hyfp_xmax, hxfp_nbins, hxfp_xmin, hxfp_xmax);
  H_exfp_vs_eyfp  = new TH2F("H_exfp_vs_eyfp", Form("%s  X_{fp} vs. Y_{fp}; Y_{fp} [cm]; X_{fp} [cm]", e_arm_name.Data()),  eyfp_nbins, eyfp_xmin, eyfp_xmax, exfp_nbins, exfp_xmin, exfp_xmax);

  
  //Add ACCP Histos to TList
  accp_HList->Add( H_exfp       );
  accp_HList->Add( H_eyfp       );
  accp_HList->Add( H_expfp      );
  accp_HList->Add( H_eypfp      );

  accp_HList->Add( H_eytar       );
  accp_HList->Add( H_exptar      );
  accp_HList->Add( H_eyptar      );
  accp_HList->Add( H_edelta      );
  
  accp_HList->Add( H_hxfp       );
  accp_HList->Add( H_hyfp       );
  accp_HList->Add( H_hxpfp      );
  accp_HList->Add( H_hypfp      );
  
  accp_HList->Add( H_hytar       );
  accp_HList->Add( H_hxptar      );
  accp_HList->Add( H_hyptar      );
  accp_HList->Add( H_hdelta      );

  accp_HList->Add( H_htar_x       );
  accp_HList->Add( H_htar_y       );
  accp_HList->Add( H_htar_z       );
  accp_HList->Add( H_etar_x       );
  accp_HList->Add( H_etar_y       );
  accp_HList->Add( H_etar_z       );
  accp_HList->Add( H_ztar_diff    );

  accp_HList->Add( H_hXColl      );
  accp_HList->Add( H_hYColl      );
  accp_HList->Add( H_eXColl      );
  accp_HList->Add( H_eYColl      );

  accp_HList->Add( H_hXColl_vs_hYColl  );
  accp_HList->Add( H_eXColl_vs_eYColl  );
  
  accp_HList->Add( H_hxfp_vs_hyfp  );
  accp_HList->Add( H_exfp_vs_eyfp  );

  
}
//_______________________________________________________________________________
void baseAnalyzer::ReadScalerTree()
{
  cout << "Calling Base ReadScalerTree()  " << endl;
  cout << Form("Using %s ", bcm_type.Data()) << endl;

  //Read ROOTfile
  inROOT = new TFile(data_InputFileName.Data(), "READ");
  
  
  //Get the SHMS scaler tree (HMS Scaler tree should be identical copy in coin mode. It is just used for cross check)
  scaler_tree = (TTree*)inROOT->Get(scl_tree_name.Data());
  scal_entries = scaler_tree->GetEntries();
  evt_flag_bcm = new Int_t[scal_entries]; //store 0 or 1, to determine which scaler read passed cut
  scal_evt_num = new Int_t[scal_entries]; //store event associated with scaler read
  
  if(daq_mode=="coin")
    {
      //SetBranchAddress
      scaler_tree->SetBranchAddress("evNumber", &Scal_evNum);
      scaler_tree->SetBranchAddress(Form("P.%s.scalerCharge", bcm_type.Data()), &Scal_BCM_charge); 
      scaler_tree->SetBranchAddress(Form("P.%s.scalerCurrent", bcm_type.Data()), &Scal_BCM_current); 
      scaler_tree->SetBranchAddress("P.1MHz.scalerTime",&Scal_time);
      scaler_tree->SetBranchAddress("P.S1X.scaler",&S1X_scaler);  
      scaler_tree->SetBranchAddress("P.pTRIG1.scaler",&TRIG1_scaler);
      scaler_tree->SetBranchAddress("P.pTRIG2.scaler",&TRIG2_scaler);
      scaler_tree->SetBranchAddress("P.pTRIG3.scaler",&TRIG3_scaler);
      scaler_tree->SetBranchAddress("P.pTRIG4.scaler",&TRIG4_scaler);
      scaler_tree->SetBranchAddress("P.pTRIG5.scaler",&TRIG5_scaler);
      scaler_tree->SetBranchAddress("P.pTRIG6.scaler",&TRIG6_scaler);
      scaler_tree->SetBranchAddress("P.EDTM.scaler",  &EDTM_scaler);
    }
  else if(daq_mode=="singles")
    {
      //SetBranchAddress
      scaler_tree->SetBranchAddress("evNumber", &Scal_evNum);
      scaler_tree->SetBranchAddress(Form("%s.%s.scalerCharge", eArm.Data(), bcm_type.Data()),  &Scal_BCM_charge); 
      scaler_tree->SetBranchAddress(Form("%s.%s.scalerCurrent", eArm.Data(), bcm_type.Data()), &Scal_BCM_current); 
      scaler_tree->SetBranchAddress(Form("%s.1MHz.scalerTime", eArm.Data()),                   &Scal_time);
      scaler_tree->SetBranchAddress(Form("%s.S1X.scaler", eArm.Data()),                        &S1X_scaler);  
      scaler_tree->SetBranchAddress(Form("%s.%sTRIG1.scaler",eArm.Data(), e_arm.Data() ),      &TRIG1_scaler);
      scaler_tree->SetBranchAddress(Form("%s.%sTRIG2.scaler",eArm.Data(), e_arm.Data() ),      &TRIG2_scaler);
      scaler_tree->SetBranchAddress(Form("%s.%sTRIG3.scaler",eArm.Data(), e_arm.Data() ),      &TRIG3_scaler);
      scaler_tree->SetBranchAddress(Form("%s.%sTRIG4.scaler",eArm.Data(), e_arm.Data() ),      &TRIG4_scaler);
      scaler_tree->SetBranchAddress(Form("%s.%sTRIG5.scaler",eArm.Data(), e_arm.Data() ),      &TRIG5_scaler);
      scaler_tree->SetBranchAddress(Form("%s.%sTRIG6.scaler",eArm.Data(), e_arm.Data() ),      &TRIG6_scaler);
      scaler_tree->SetBranchAddress(Form("%s.EDTM.scaler",eArm.Data() ),                       &EDTM_scaler);
    }
  
   
  
  cout << "Ending ReadScalerTree() . . . " << endl;
  
}

//_______________________________________________________________________________
void baseAnalyzer::ScalerEventLoop()
{
  
  cout << "Calling Base ScalerEventLoop() " << endl;
  

  //Scaler reads loop. ith scaler read
  for (int i = 0; i < scal_entries; i++) 
    {
      /*(NOTE: Each scaler read is associated with as specific event number
	as (scaler read 1-> event 1000,  scaler read 2 -> event 2300, ...)
	This means events up to 1000 correspond to scaler read 1, ...*/
      
      scaler_tree->GetEntry(i);
      
      //Set Event Flag to FALSE (default)
      evt_flag_bcm[i] = 0;
      
      //Store event associated with scaler read
      scal_evt_num[i] = Scal_evNum;
      
      //Store Cumulative Quantities
      total_charge_bcm = Scal_BCM_charge;
      total_time = Scal_time;
      total_s1x_scaler = S1X_scaler;
      total_trig1_scaler = TRIG1_scaler;
      total_trig2_scaler = TRIG2_scaler;
      total_trig3_scaler = TRIG3_scaler;
      total_trig4_scaler = TRIG4_scaler;
      total_trig5_scaler = TRIG5_scaler;
      total_trig6_scaler = TRIG6_scaler;
      total_edtm_scaler = EDTM_scaler;

    
      //Check If BCM Beam Current in Between Reads is Over Threshold
      // if(abs(Scal_BCM_current-set_current)<bcm_thrs)
      if(Scal_BCM_current>bcm_thrs)
	{
	  
	  //Turn Event Flag ON, if beam current is within threshold
	  evt_flag_bcm[i] = 1;
	  
	  //Store Quantities that Passed the Current Threshold
	  total_time_bcm_cut = total_time_bcm_cut + (Scal_time - prev_time);
	  total_charge_bcm_cut = total_charge_bcm_cut + (Scal_BCM_charge - prev_charge_bcm);  
	  total_s1x_scaler_bcm_cut = total_s1x_scaler_bcm_cut + (S1X_scaler-prev_s1x_scaler);
	  total_trig1_scaler_bcm_cut = total_trig1_scaler_bcm_cut + (TRIG1_scaler-prev_trig1_scaler);
	  total_trig2_scaler_bcm_cut = total_trig2_scaler_bcm_cut + (TRIG2_scaler-prev_trig2_scaler);
	  total_trig3_scaler_bcm_cut = total_trig3_scaler_bcm_cut + (TRIG3_scaler-prev_trig3_scaler);
	  total_trig4_scaler_bcm_cut = total_trig4_scaler_bcm_cut + (TRIG4_scaler-prev_trig4_scaler);
	  total_trig5_scaler_bcm_cut = total_trig5_scaler_bcm_cut + (TRIG5_scaler-prev_trig5_scaler);
	  total_trig6_scaler_bcm_cut = total_trig6_scaler_bcm_cut + (TRIG6_scaler-prev_trig6_scaler);
	  total_edtm_scaler_bcm_cut = total_edtm_scaler_bcm_cut + (EDTM_scaler - prev_edtm_scaler);

	} //End BCM Current Cut

      //Previous Scaler Reads (Necessary to Take Average between S-1 and S scaler reads, to get values in between)
      prev_time = Scal_time;
      prev_charge_bcm = Scal_BCM_charge;
      prev_s1x_scaler = S1X_scaler;
      prev_trig1_scaler = TRIG1_scaler;
      prev_trig2_scaler = TRIG2_scaler;
      prev_trig3_scaler = TRIG3_scaler;
      prev_trig4_scaler = TRIG4_scaler;
      prev_trig5_scaler = TRIG5_scaler;
      prev_trig6_scaler = TRIG6_scaler;
      prev_edtm_scaler = EDTM_scaler;
      
      // print every 100000 scaler reads
      if (i % 100000 == 0){  
	cout << "ScalerEventLoop(PASS2): " << std::setprecision(2) << double(i) / scal_entries * 100. << "  % " << std::flush << "\r";
      }
    
    } //End Scaler Read Loop
   
  //Subtract EDTM counts from trigger scalers
  total_s1x_scaler_bcm_cut = total_s1x_scaler_bcm_cut - total_edtm_scaler_bcm_cut;
  total_trig1_scaler_bcm_cut = total_trig1_scaler_bcm_cut - total_edtm_scaler_bcm_cut;
  total_trig2_scaler_bcm_cut = total_trig2_scaler_bcm_cut - total_edtm_scaler_bcm_cut;
  total_trig3_scaler_bcm_cut = total_trig3_scaler_bcm_cut - total_edtm_scaler_bcm_cut;
  total_trig4_scaler_bcm_cut = total_trig4_scaler_bcm_cut - total_edtm_scaler_bcm_cut;
  total_trig5_scaler_bcm_cut = total_trig5_scaler_bcm_cut - total_edtm_scaler_bcm_cut;
  total_trig6_scaler_bcm_cut = total_trig6_scaler_bcm_cut - total_edtm_scaler_bcm_cut;

  //Calculate Trigger Rates (EDTM subtracted already)
  S1XscalerRate_bcm_cut = total_s1x_scaler_bcm_cut / total_time_bcm_cut;
  TRIG1scalerRate_bcm_cut = total_trig1_scaler_bcm_cut / total_time_bcm_cut;
  TRIG2scalerRate_bcm_cut = total_trig2_scaler_bcm_cut / total_time_bcm_cut;
  TRIG3scalerRate_bcm_cut = total_trig3_scaler_bcm_cut / total_time_bcm_cut;
  TRIG4scalerRate_bcm_cut = total_trig4_scaler_bcm_cut / total_time_bcm_cut;
  TRIG5scalerRate_bcm_cut = total_trig5_scaler_bcm_cut / total_time_bcm_cut;
  TRIG6scalerRate_bcm_cut = total_trig6_scaler_bcm_cut / total_time_bcm_cut;
  EDTMscalerRate_bcm_cut =  total_edtm_scaler_bcm_cut / total_time_bcm_cut;
    
  cout << "Ending ScalerEventLoop() . . . " << endl;
  
}

//_______________________________________________________________________________
void baseAnalyzer::ReadTree()
{
  cout << "Calling Base ReadTree()  " << endl;

  
  if(analysis=="data")
    {
      
      cout << "Analyzing DATA . . . " << endl;

      //Read ROOTfile
      inROOT = new TFile(data_InputFileName.Data(), "READ");
      
      //Get the data tree
      tree = (TTree*)inROOT->Get("T");
      nentries = tree->GetEntries();
            
      //---------------SetBranchAddress-----------------

      // Global Variables
      tree->SetBranchAddress("g.evtyp",&gevtyp);
      tree->SetBranchAddress("g.evnum",&gevnum);

      //This is meant to be used by the derived class, helicityAnalyzer.cpp (see 'helicity_flag' in main_controls.cpp)
      if(helicity_flag)
	{	  
	  tree->SetBranchAddress("T.helicity.cycle",         &hel_cycle);  //Helicity Cycle
	  tree->SetBranchAddress("T.helicity.hel",           &hel);        //actual helicity for event
	  tree->SetBranchAddress("T.helicity.helpred",       &hel_pred);   //predicted reported helicity for event
	  tree->SetBranchAddress("T.helicity.helrep",        &hel_rep);    //reported helicity for event
	  tree->SetBranchAddress("T.helicity.mps",           &hel_mps);    //In MPS blanking period (helicity unknown or undetermined)
	  tree->SetBranchAddress("T.helicity.nqrt",          &hel_nqrt);   //position of helicity in quartet
	  tree->SetBranchAddress("T.helicity.pcheck",        &hel_pcheck);  //Period check
	  tree->SetBranchAddress("T.helicity.qrt",           &hel_qrt);     //Last cycle of quartet
	  
	}
      
      if(daq_mode=="coin"){

	//Coincidence Time
	tree->SetBranchAddress("CTime.epCoinTime_ROC2",  &epCoinTime);
	tree->SetBranchAddress("CTime.eKCoinTime_ROC2",  &eKCoinTime);
	tree->SetBranchAddress("CTime.ePiCoinTime_ROC2", &ePiCoinTime);
	
	// Trigger Detector 
	tree->SetBranchAddress("T.coin.pTRIG1_ROC2_tdcTimeRaw",&TRIG1_tdcTimeRaw);
	tree->SetBranchAddress("T.coin.pTRIG2_ROC2_tdcTimeRaw",&TRIG2_tdcTimeRaw);
	tree->SetBranchAddress("T.coin.pTRIG3_ROC2_tdcTimeRaw",&TRIG3_tdcTimeRaw);
	tree->SetBranchAddress("T.coin.pTRIG4_ROC2_tdcTimeRaw",&TRIG4_tdcTimeRaw);
	tree->SetBranchAddress("T.coin.pTRIG5_ROC2_tdcTimeRaw",&TRIG5_tdcTimeRaw);
	tree->SetBranchAddress("T.coin.pTRIG6_ROC2_tdcTimeRaw",&TRIG6_tdcTimeRaw);
	tree->SetBranchAddress("T.coin.pEDTM_tdcTimeRaw",&EDTM_tdcTimeRaw);

	//------------------------------------
	//-----Kinematics Leaf Variables------
	//------------------------------------
	//Primary Kinematics (electron kinematics)
	tree->SetBranchAddress(Form("%s.kin.primary.scat_ang_rad", eArm.Data()),&th_e);
	tree->SetBranchAddress(Form("%s.kin.primary.W", eArm.Data()),&W);
	tree->SetBranchAddress(Form("%s.kin.primary.W2", eArm.Data()),&W2);
	tree->SetBranchAddress(Form("%s.kin.primary.Q2", eArm.Data()),&Q2);
	tree->SetBranchAddress(Form("%s.kin.primary.x_bj", eArm.Data()),&X);
	tree->SetBranchAddress(Form("%s.kin.primary.nu", eArm.Data()),&nu);
	tree->SetBranchAddress(Form("%s.kin.primary.q3m", eArm.Data()),&q);
	tree->SetBranchAddress(Form("%s.kin.primary.q_x", eArm.Data()),&qx);
	tree->SetBranchAddress(Form("%s.kin.primary.q_y", eArm.Data()),&qy);
	tree->SetBranchAddress(Form("%s.kin.primary.q_z", eArm.Data()),&qz);
	tree->SetBranchAddress(Form("%s.kin.primary.th_q", eArm.Data()),&th_q);
	tree->SetBranchAddress(Form("%s.kin.primary.ph_q", eArm.Data()),&ph_q);
	tree->SetBranchAddress(Form("%s.kin.primary.epsilon", eArm.Data()),&epsilon);
	
	//Secondary Kinematics (hadron kinematics)
	tree->SetBranchAddress(Form("%s.kin.secondary.emiss", hArm.Data()),&Em);
	tree->SetBranchAddress(Form("%s.kin.secondary.emiss_nuc", hArm.Data()),&Em_nuc);     
	tree->SetBranchAddress(Form("%s.kin.secondary.pmiss", hArm.Data()),&Pm);
	tree->SetBranchAddress(Form("%s.kin.secondary.Prec_x", hArm.Data()),&Pmx_lab);       
	tree->SetBranchAddress(Form("%s.kin.secondary.Prec_y", hArm.Data()),&Pmy_lab);   
	tree->SetBranchAddress(Form("%s.kin.secondary.Prec_z", hArm.Data()),&Pmz_lab);   
	tree->SetBranchAddress(Form("%s.kin.secondary.pmiss_x", hArm.Data()),&Pmx_q);        
	tree->SetBranchAddress(Form("%s.kin.secondary.pmiss_y", hArm.Data()),&Pmy_q);   
	tree->SetBranchAddress(Form("%s.kin.secondary.pmiss_z", hArm.Data()),&Pmz_q);   
	tree->SetBranchAddress(Form("%s.kin.secondary.tx", hArm.Data()),&Tx);                 
	tree->SetBranchAddress(Form("%s.kin.secondary.tb", hArm.Data()),&Tr);                 
	tree->SetBranchAddress(Form("%s.kin.secondary.Mrecoil", hArm.Data()),&MM);      
	tree->SetBranchAddress(Form("%s.kin.secondary.th_xq", hArm.Data()),&th_xq);            
	tree->SetBranchAddress(Form("%s.kin.secondary.th_bq", hArm.Data()),&th_rq);            
	tree->SetBranchAddress(Form("%s.kin.secondary.ph_xq", hArm.Data()),&ph_xq);                           
	tree->SetBranchAddress(Form("%s.kin.secondary.ph_bq", hArm.Data()),&ph_rq);                            
	tree->SetBranchAddress(Form("%s.kin.secondary.xangle", hArm.Data()),&xangle);    
	tree->SetBranchAddress(Form("%s.kin.secondary.tx_cm", hArm.Data()),&Tx_cm);                 
	tree->SetBranchAddress(Form("%s.kin.secondary.tb_cm", hArm.Data()),&Tr_cm);
	tree->SetBranchAddress(Form("%s.kin.secondary.thx_cm", hArm.Data()),&th_xq_cm);            
	tree->SetBranchAddress(Form("%s.kin.secondary.thb_cm", hArm.Data()),&th_rq_cm);            
	tree->SetBranchAddress(Form("%s.kin.secondary.phx_cm", hArm.Data()),&ph_xq_cm);                           
	tree->SetBranchAddress(Form("%s.kin.secondary.phb_cm", hArm.Data()),&ph_rq_cm);
	tree->SetBranchAddress(Form("%s.kin.secondary.t_tot_cm", hArm.Data()),&Ttot_cm);
	tree->SetBranchAddress(Form("%s.kin.secondary.MandelS", hArm.Data()),&MandelS);
	tree->SetBranchAddress(Form("%s.kin.secondary.MandelT", hArm.Data()),&MandelT);
	tree->SetBranchAddress(Form("%s.kin.secondary.MandelU", hArm.Data()),&MandelU);

	//Additional Kinematics (not in primary/secondary modules)
	tree->SetBranchAddress(Form("%s.gtr.p",  hArm.Data()), &Pf);
	
	//----Hadron Arm Focal Plane----- 
	tree->SetBranchAddress(Form("%s.dc.x_fp",  hArm.Data()), &h_xfp);
	tree->SetBranchAddress(Form("%s.dc.xp_fp", hArm.Data()), &h_xpfp);
	tree->SetBranchAddress(Form("%s.dc.y_fp",  hArm.Data()), &h_yfp);
	tree->SetBranchAddress(Form("%s.dc.yp_fp", hArm.Data()), &h_ypfp);
	
	//----Hadron Arm Reconstructed-----
	tree->SetBranchAddress(Form("%s.gtr.y",  hArm.Data()), &h_ytar);
	tree->SetBranchAddress(Form("%s.gtr.ph", hArm.Data()), &h_yptar);
	tree->SetBranchAddress(Form("%s.gtr.th", hArm.Data()), &h_xptar);
	tree->SetBranchAddress(Form("%s.gtr.dp", hArm.Data()), &h_delta);

	//----Target Quantities----
	//(tarx, tary, tarz) in Hall Coord. System
	tree->SetBranchAddress(Form("%s.react.x", hArm.Data()), &htar_x);
	tree->SetBranchAddress(Form("%s.react.y", hArm.Data()), &htar_y);
	tree->SetBranchAddress(Form("%s.react.z", hArm.Data()), &htar_z);

	//----Collimator Quantities-----
	tree->SetBranchAddress(Form("%s.extcor.xsieve", hArm.Data()),&hXColl);
	tree->SetBranchAddress(Form("%s.extcor.ysieve", hArm.Data()),&hYColl);
	
      }
      
      if(daq_mode=="singles"){
	// Trigger Detector 
	tree->SetBranchAddress(Form("T.%s.%sTRIG1_tdcTimeRaw", daq.Data(), e_arm.Data() ),&TRIG1_tdcTimeRaw);
	tree->SetBranchAddress(Form("T.%s.%sTRIG2_tdcTimeRaw", daq.Data(), e_arm.Data() ),&TRIG2_tdcTimeRaw);
	tree->SetBranchAddress(Form("T.%s.%sTRIG3_tdcTimeRaw", daq.Data(), e_arm.Data() ),&TRIG3_tdcTimeRaw);
	tree->SetBranchAddress(Form("T.%s.%sTRIG4_tdcTimeRaw", daq.Data(), e_arm.Data() ),&TRIG4_tdcTimeRaw);
	tree->SetBranchAddress(Form("T.%s.%sTRIG5_tdcTimeRaw", daq.Data(), e_arm.Data() ),&TRIG5_tdcTimeRaw);
	tree->SetBranchAddress(Form("T.%s.%sTRIG6_tdcTimeRaw", daq.Data(), e_arm.Data() ),&TRIG6_tdcTimeRaw);
	tree->SetBranchAddress(Form("T.%s.%sEDTM_tdcTimeRaw", daq.Data(), e_arm.Data()),&EDTM_tdcTimeRaw);

	//------------------------------------
	//-----Kinematics Leaf Variables------
	//------------------------------------
	//Primary Kinematics (electron kinematics) on single-arm mode
	tree->SetBranchAddress(Form("%s.kin.scat_ang_rad", eArm.Data()),&th_e);
	tree->SetBranchAddress(Form("%s.kin.W", eArm.Data()),&W);
	tree->SetBranchAddress(Form("%s.kin.W2", eArm.Data()),&W2);
	tree->SetBranchAddress(Form("%s.kin.Q2", eArm.Data()),&Q2);
	tree->SetBranchAddress(Form("%s.kin.x_bj", eArm.Data()),&X);
	tree->SetBranchAddress(Form("%s.kin.nu", eArm.Data()),&nu);
	tree->SetBranchAddress(Form("%s.kin.q3m", eArm.Data()),&q);
	tree->SetBranchAddress(Form("%s.kin.q_x", eArm.Data()),&qx);
	tree->SetBranchAddress(Form("%s.kin.q_y", eArm.Data()),&qy);
	tree->SetBranchAddress(Form("%s.kin.q_z", eArm.Data()),&qz);
	tree->SetBranchAddress(Form("%s.kin.th_q", eArm.Data()),&th_q);
	tree->SetBranchAddress(Form("%s.kin.ph_q", eArm.Data()),&ph_q);
	tree->SetBranchAddress(Form("%s.kin.epsilon", eArm.Data()),&epsilon);
	
      }
      

      if(daq_mode=="coin" || (daq_mode=="singles" && e_arm_name=="HMS") )
	{
	  //-------------------------------
	  //------PID Leaf Variables-------
	  //-------------------------------
	  
	  //HMS DETECTORS
	  tree->SetBranchAddress("H.cer.npeSum",         &hcer_npesum);
	  tree->SetBranchAddress("H.cal.etotnorm",       &hcal_etotnorm);
	  tree->SetBranchAddress("H.cal.etottracknorm",  &hcal_etottracknorm);
	  tree->SetBranchAddress("H.hod.betanotrack",    &hhod_beta_ntrk);
	  tree->SetBranchAddress("H.hod.beta",           &hhod_beta);
	  tree->SetBranchAddress("H.gtr.beta",           &hhod_gtr_beta);
	  tree->SetBranchAddress("H.hod.goodscinhit",    &hhod_GoodScinHit);    
	  tree->SetBranchAddress("H.dc.ntrack",          &hdc_ntrack);

	  //-----------------------------------
	  //-----Acceptance Leaf Variables-----
	  //-----------------------------------

	  //Additional Kinematics (not in primary/secondary modules)
	  tree->SetBranchAddress(Form("%s.gtr.p",  eArm.Data()), &kf);
	  
	  //----Electron Arm Focal Plane---- 
	  tree->SetBranchAddress(Form("%s.dc.x_fp",  eArm.Data()), &e_xfp);
	  tree->SetBranchAddress(Form("%s.dc.xp_fp", eArm.Data()), &e_xpfp);
	  tree->SetBranchAddress(Form("%s.dc.y_fp",  eArm.Data()), &e_yfp);
	  tree->SetBranchAddress(Form("%s.dc.yp_fp", eArm.Data()), &e_ypfp);
	  
	  //----Electron Arm Reconstructed----- 
	  tree->SetBranchAddress(Form("%s.gtr.y",  eArm.Data()), &e_ytar);
	  tree->SetBranchAddress(Form("%s.gtr.ph", eArm.Data()), &e_yptar);
	  tree->SetBranchAddress(Form("%s.gtr.th", eArm.Data()), &e_xptar);
	  tree->SetBranchAddress(Form("%s.gtr.dp", eArm.Data()), &e_delta);
	  
	  //----Target Quantities----
	  //(tarx, tary, tarz) in Hall Coord. System      
	  tree->SetBranchAddress(Form("%s.react.x", eArm.Data()), &etar_x);
	  tree->SetBranchAddress(Form("%s.react.y", eArm.Data()), &etar_y);
	  tree->SetBranchAddress(Form("%s.react.z", eArm.Data()), &etar_z);

	  //----Collimator Quantities-----
	  tree->SetBranchAddress(Form("%s.extcor.xsieve", eArm.Data()),&eXColl);
	  tree->SetBranchAddress(Form("%s.extcor.ysieve", eArm.Data()),&eYColl);
	  
	}
      
      if(daq_mode=="coin" || (daq_mode=="singles" && e_arm_name=="SHMS") )
	{
	  //-------------------------------
	  //------PID Leaf Variables-------
	  //-------------------------------

	  //SHMS DETECTORS
	  tree->SetBranchAddress("P.ngcer.npeSum",       &pngcer_npesum);
	  tree->SetBranchAddress("P.hgcer.npeSum",       &phgcer_npesum);
	  tree->SetBranchAddress("P.aero.npeSum",        &paero_npesum);
	  tree->SetBranchAddress("P.cal.etotnorm",       &pcal_etotnorm);
	  tree->SetBranchAddress("P.cal.etottracknorm",  &pcal_etottracknorm);
	  tree->SetBranchAddress("P.hod.betanotrack",    &phod_beta_ntrk);
	  tree->SetBranchAddress("P.hod.beta",           &phod_beta);
	  tree->SetBranchAddress("P.gtr.beta",           &phod_gtr_beta);
	  tree->SetBranchAddress("P.hod.goodscinhit",    &phod_GoodScinHit);    
	  tree->SetBranchAddress("P.dc.ntrack",          &pdc_ntrack);

	  //-----------------------------------
	  //-----Acceptance Leaf Variables-----
	  //-----------------------------------

	  //Additional Kinematics (not in primary/secondary modules)
	  tree->SetBranchAddress(Form("%s.gtr.p",  eArm.Data()), &kf);
	  
	  //----Electron Arm Focal Plane---- 
	  tree->SetBranchAddress(Form("%s.dc.x_fp",  eArm.Data()), &e_xfp);
	  tree->SetBranchAddress(Form("%s.dc.xp_fp", eArm.Data()), &e_xpfp);
	  tree->SetBranchAddress(Form("%s.dc.y_fp",  eArm.Data()), &e_yfp);
	  tree->SetBranchAddress(Form("%s.dc.yp_fp", eArm.Data()), &e_ypfp);
	  
	  //----Electron Arm Reconstructed----- 
	  tree->SetBranchAddress(Form("%s.gtr.y",  eArm.Data()), &e_ytar);
	  tree->SetBranchAddress(Form("%s.gtr.ph", eArm.Data()), &e_yptar);
	  tree->SetBranchAddress(Form("%s.gtr.th", eArm.Data()), &e_xptar);
	  tree->SetBranchAddress(Form("%s.gtr.dp", eArm.Data()), &e_delta);

	  //----Target Quantities----
	  //(tarx, tary, tarz) in Hall Coord. System      
	  tree->SetBranchAddress(Form("%s.react.x", eArm.Data()), &etar_x);
	  tree->SetBranchAddress(Form("%s.react.y", eArm.Data()), &etar_y);
	  tree->SetBranchAddress(Form("%s.react.z", eArm.Data()), &etar_z);

	  //----Collimator Quantities-----
	  tree->SetBranchAddress(Form("%s.extcor.xsieve", eArm.Data()),&eXColl);
	  tree->SetBranchAddress(Form("%s.extcor.ysieve", eArm.Data()),&eYColl);
	  
	}


    } //END DATA SET BRANCH ADDRESS

  else if(analysis=="simc")
    {
      cout << "SIMC ANALYSIS C++ CODE HAS NOT BEEN DONE YET ! ! !" << endl;
    } //END SIMC SET BRANCH ADDRESS


}

//_______________________________________________________________________________
void baseAnalyzer::EventLoop()
{

  cout << "Calling Base EventLoop() . . . " << endl;
  
  //Loop over Events
  
  if(analysis=="data")
    {
      cout << "Loop over Data Events | nentries -->  " << nentries << endl;

      for(int ientry=0; ientry<nentries; ientry++)
	{
	  
	  tree->GetEntry(ientry);

     
	  //--------------CALCULATED KINEMATICS VARIABLES (IF THEY ARE NOT ALREADY DONE IN HCANA)-----------

	  th_x = xangle - th_e;  //hadron arm central angle for each particle
	  MM2 = MM*MM;           //Missing Mass Squared
 	  ztar_diff = htar_z - etar_z;  //reaction vertex z difference
	  
	  
	  
	  //--------------DEFINE CUTS--------------------
	  
	  //CUTS USED IN EDTM LIVE TIME CALCULATION
	  c_noedtm = EDTM_tdcTimeRaw == 0.;
	  c_edtm   = EDTM_tdcTimeRaw  > 0.;
	  c_trig1  = TRIG1_tdcTimeRaw > 0.;
	  c_trig2  = TRIG2_tdcTimeRaw > 0.;
	  c_trig3  = TRIG3_tdcTimeRaw > 0.;
	  c_trig4  = TRIG4_tdcTimeRaw > 0.;
	  c_trig5  = TRIG5_tdcTimeRaw > 0.;
	  c_trig6  = TRIG6_tdcTimeRaw > 0.;

	  //=====CUTS USED IN TRACKING EFFICIENCY CALCULATION=====

	  //CUTS: HMS TRACKING EFFICIENCY (May be for e- or hadrons, depending on the limits set in the input file)

	  //Require at least a minimum number of track(s) 
	  if(hdc_ntrk_cut_flag){c_hdc_ntrk = hdc_ntrack >= c_hdc_ntrk_min;}
	  else{
	    cout <<
	      "********************************\n"
	      "TRACKING EFFICIENCY ERROR: \n"
	      "Must set hdc_ntrk_cut_flag = 1 \n"
	      "See " <<  Form("%s", input_CutFileName.Data()) << "\n"
	      "********************************"<< endl;
	    
	    //"*********************************************************" << endl;
	    gSystem->Exit(0);}

	  //Require a "Good Scintillator Hit" in the Fiducial Hodoscopoe Region
	  if(hScinGood_cut_flag){c_hScinGood = hhod_GoodScinHit==1;}
	  else{c_hScinGood=1;} //1 means allow all events (i.e., do not put hScinGood cut other than 1)

	  //Require HMS Cherenkov Cut (for electron or hadron selection)
	  if(hcer_cut_flag){c_hcer_NPE_Sum = hcer_npesum >= c_hnpeSum_min && hcer_npesum <= c_hnpeSum_max;}
	  else{c_hcer_NPE_Sum=1;}

	  //Require HMS Calorimeter Cut (for additional electron or hadron selection)
	  if(hetotnorm_cut_flag){c_hetotnorm = hcal_etotnorm >= c_hetotnorm_min && hcal_etotnorm <= c_hetotnorm_max;}
	  else{c_hetotnorm=1;}

	  //Require HMS Hodoscope Beta Cut (no track-biased) (for electron or hadron selection)
	  if(hBeta_notrk_cut_flag){c_hBeta_notrk = hhod_beta_ntrk >= c_hBetaNtrk_min && hhod_beta_ntrk <= c_hBetaNtrk_max;}
	  else{c_hBeta_notrk=1;}

	  //electrons (or hadrons) that 'SHOULD' have passed the cuts to form a track
	  good_hms_should = c_hScinGood && c_hcer_NPE_Sum && c_hetotnorm && c_hBeta_notrk;

	  //electrons (or hadrons) that 'DID' passed the cuts to form a track
	  good_hms_did = c_hdc_ntrk && good_hms_should;


	  //CUTS: SHMS TRACKING EFFICIENCY (May be for e- or hadrons, depending on the limits set in the input file)

	  //Require at least a minimum number of track(s) 
	  if(pdc_ntrk_cut_flag){c_pdc_ntrk = pdc_ntrack >= c_pdc_ntrk_min;}
	  else{
	    cout <<
	      "********************************\n"
	      "TRACKING EFFICIENCY ERROR: \n"
	      "Must set pdc_ntrk_cut_flag = 1 \n"
	      "See " <<  Form("%s", input_CutFileName.Data()) << "\n"
	      "********************************"<< endl;
	    
	    //"*********************************************************" << endl;
	    gSystem->Exit(0);}

	  //Require a "Good Scintillator Hit" in the Fiducial Hodoscopoe Region
	  if(pScinGood_cut_flag){c_pScinGood = phod_GoodScinHit==1;}
	  else{c_pScinGood=1;} //1 means allow all events (do not put hScinGood cut)

	  //Require SHMS Noble Gas Cherenkov Cut (for electron or hadron selection)
	  if(pngcer_cut_flag){c_pngcer_NPE_Sum = pngcer_npesum >= c_pngcer_npeSum_min &&  pngcer_npesum <= c_pngcer_npeSum_max;}
	  else{c_pngcer_NPE_Sum=1;}

	  //Require SHMS Heavy Gas Cherenkov Cut (for electron or hadron selection)
	  if(phgcer_cut_flag){c_phgcer_NPE_Sum = phgcer_npesum >= c_phgcer_npeSum_min &&  phgcer_npesum <= c_phgcer_npeSum_max;}
	  else{c_phgcer_NPE_Sum=1;}

	  //WHAT ABOUT AN AEROGEL CUT IN THE TRACKING EFFICIENCY?
	  
	  //Require SHMS Calorimeter Cut (for additional electron or hadron selection)
	  if(petotnorm_cut_flag){c_petotnorm = pcal_etotnorm >= c_petotnorm_min && pcal_etotnorm <= c_petotnorm_max;}
	  else{c_petotnorm=1;}

	  //Require SHMS Hodoscope Beta Cut (no track-biased) (for electron or hadron selection)
	  if(pBeta_notrk_cut_flag){c_pBeta_notrk = phod_beta_ntrk >= c_pBetaNtrk_min && phod_beta_ntrk <= c_pBetaNtrk_max;}
	  else{c_pBeta_notrk=1;}

	  //electrons (or hadrons) that 'SHOULD' have passed the cuts to form a track
	  good_shms_should = c_pScinGood && c_pngcer_NPE_Sum && c_phgcer_NPE_Sum && c_petotnorm && c_pBeta_notrk;

	  //electrons (or hadrons) that 'DID' passed the cuts to form a track
	  good_shms_did = c_pdc_ntrk && good_shms_should;

	  //=====END: CUTS USED IN TRACKING EFFICIENCY CALCULATION=====
	  
	  
	  //====DATA ANALYSIS CUTS (MUST BE EXACTLY SAME AS SIMC, except PID CUTS on detectors)====

	  //----PID Cuts---- (SPECIFIC TO DATA)

	  // kaon coincidence time cut
	  if(eKctime_pidCut_flag) {cpid_eK_ctime = eKCoinTime>=cpid_eKctime_min && eKCoinTime<=cpid_eKctime_max;}
	  else{cpid_eK_ctime=1;}

	  // pion coincidence time cut
	  if(ePictime_pidCut_flag) {cpid_ePi_ctime = ePiCoinTime>=cpid_ePictime_min && ePiCoinTime<=cpid_ePictime_max;}
	  else{cpid_ePi_ctime=1;}

	  // proton coincidence time cut
	  if(ePctime_pidCut_flag) {cpid_eP_ctime = epCoinTime>=cpid_ePctime_min && epCoinTime<=cpid_ePctime_max;}
	  else{cpid_eP_ctime=1;}

	  //SHMS calorimeter total normalized track energy
	  if(petot_trkNorm_pidCut_flag) {cpid_petot_trkNorm = pcal_etottracknorm>=cpid_petot_trkNorm_min && pcal_etottracknorm<=cpid_petot_trkNorm_max;}
	  else{cpid_petot_trkNorm=1;}

	  //SHMS Noble Gas Cherenkov
	  if(pngcer_pidCut_flag) {cpid_pngcer_NPE_Sum = pngcer_npesum>=cpid_pngcer_npeSum_min && pngcer_npesum<=cpid_pngcer_npeSum_max;}
	  else{cpid_pngcer_NPE_Sum=1;}

	  //SHMS Heavy Gas Cherenkov
	  if(phgcer_pidCut_flag) {cpid_phgcer_NPE_Sum = phgcer_npesum>=cpid_phgcer_npeSum_min && phgcer_npesum<=cpid_phgcer_npeSum_max;}
	  else{cpid_phgcer_NPE_Sum=1;}

	  //SHMS Aerogel Cherenkov
	  if(paero_pidCut_flag) {cpid_paero_NPE_Sum = paero_npesum>=cpid_paero_npeSum_min && paero_npesum<=cpid_paero_npeSum_max;}
	  else{cpid_paero_NPE_Sum=1;}

	  //HMS calorimeter total normalized track energy
	  if(hetot_trkNorm_pidCut_flag) {cpid_hetot_trkNorm = hcal_etottracknorm>=cpid_hetot_trkNorm_min && hcal_etottracknorm<=cpid_hetot_trkNorm_max;}
	  else{cpid_hetot_trkNorm=1;}

	  //HMS Gas Cherenkov
	  if(hcer_pidCut_flag) {cpid_hcer_NPE_Sum = hcer_npesum>=cpid_hcer_npeSum_min && hcer_npesum<=cpid_hcer_npeSum_max;}
	  else{cpid_hcer_NPE_Sum=1;}
	  
	  //----Kinematics Cuts----
	  //Q2
	  if(Q2_cut_flag){c_Q2 = Q2>=c_Q2_min && Q2<=c_Q2_max;}
	  else{c_Q2=1;}

	  //Missing Energy, Em
	  if(Em_cut_flag){c_Em = Em_nuc>=c_Em_min && Em_nuc<=c_Em_max;}
	  else{c_Em=1;}

	  //Invariant Mass, W
	  if(W_cut_flag){c_W = W>=c_W_min && W<=c_W_max;}
	  else{c_W=1;}

	  //Missing Mass, MM
	  //Kaons
	  if(MM_K_cut_flag){c_MM_K = MM>=c_MM_K_min && MM<=c_MM_K_max;}
	  else{c_MM_K=1;}
	  //Pions
	  if(MM_Pi_cut_flag){c_MM_Pi = MM>=c_MM_Pi_min && MM<=c_MM_Pi_max;}
	  else{c_MM_Pi=1;}
	  //Protons
	  if(MM_P_cut_flag){c_MM_P = MM>=c_MM_P_min && MM<=c_MM_P_max;}
	  else{c_MM_P=1;}
	  
	  //----Acceptance Cuts----
	  if(hdelta_cut_flag){c_hdelta = h_delta>=c_hdelta_min && h_delta<=c_hdelta_max;} 
	  else{c_hdelta=1;}
		  
	  if(edelta_cut_flag){c_edelta = e_delta>=c_edelta_min && e_delta<=c_edelta_max;} 
	  else{c_edelta=1;} 

	  if(ztarDiff_cut_flag){c_ztarDiff = ztar_diff>=c_ztarDiff_min && ztar_diff<=c_ztarDiff_max;} 
	  else{c_ztarDiff=1;}
	  
	  //Combine All CUTS
	  c_accpCuts = c_hdelta && c_edelta && c_ztarDiff;
	  c_pidCuts = cpid_eP_ctime && cpid_petot_trkNorm && cpid_pngcer_NPE_Sum && cpid_phgcer_NPE_Sum && cpid_paero_NPE_Sum && cpid_hetot_trkNorm && cpid_hcer_NPE_Sum;
	  c_kinCuts = c_Q2 && c_Em && c_W && c_MM_P;
	  c_baseCuts =  c_accpCuts && c_pidCuts && c_kinCuts;
	  
	  //====END: DATA ANALYSIS CUTS (MUST BE EXACTLY SAME AS SIMC)===

	  
	  //----------END DEFINE CUTS--------------------

	  //Count Accepted EDTM events (no bcm current cut : this is just to compare with counts that have bcm cuts)
	  if(c_edtm){total_edtm_accp++;}
	  
	  //Count Accepted TRIG 1-6 events (no bcm current cut : this is just to compare with counts that have bcm cuts)
	  if(c_trig1){total_trig1_accp++;}
	  if(c_trig2){total_trig2_accp++;}
	  if(c_trig3){total_trig3_accp++;}
	  if(c_trig4){total_trig4_accp++;}
	  if(c_trig5){total_trig5_accp++;}
	  if(c_trig6){total_trig6_accp++;}
	  
	  //----------------------Check If BCM Current is within limits---------------------

	  if(evt_flag_bcm[scal_read]==1)
	    {

	      //Count Accepted EDTM events (With bcm current cut: to be used in total edtm live time calculation)
	      if(c_edtm){ total_edtm_accp_bcm_cut++;}
	      
	      //Count Accepted TRIG1-6 events (without EDTM and with bcm current cut: to be used in the computer live time calculation)
	      if(c_trig1 && c_noedtm) { total_trig1_accp_bcm_cut++; }
	      if(c_trig2 && c_noedtm) { total_trig2_accp_bcm_cut++; }
	      if(c_trig3 && c_noedtm) { total_trig3_accp_bcm_cut++; }
	      if(c_trig4 && c_noedtm) { total_trig4_accp_bcm_cut++; }
	      if(c_trig5 && c_noedtm) { total_trig5_accp_bcm_cut++; }
	      if(c_trig6 && c_noedtm) { total_trig6_accp_bcm_cut++; }
	      
	      //REQUIRE "NO EDTM" CUT TO FILL DATA HISTOGRAMS
	      if(c_noedtm)
		{

		  //Calculate HMS Tracking Efficiency Components
		  if(good_hms_did){ h_did++;}
		  if(good_hms_should){ h_should++; }
		  
		  //Calculate SHMS Tracking Efficiency Components
		  if(good_shms_did){ p_did++;}
		  if(good_shms_should){ p_should++; }

		  
		  //----------------------Fill DATA Histograms-----------------------

		  if(c_baseCuts)
		    {

		      //--------------------------------------------------------------------
		      //---------HISTOGRAM CATEGORY: Particle Identification (PID)----------
		      //--------------------------------------------------------------------
		      
		      //Coincidence Time		      
		      H_ep_ctime->Fill(epCoinTime);
		      H_eK_ctime->Fill(eKCoinTime);
		      H_ePi_ctime->Fill(ePiCoinTime);

		      //Fill HMS Detectors
		      H_hCerNpeSum->Fill(hcer_npesum);
		      H_hCalEtotNorm->Fill(hcal_etotnorm);
		      H_hCalEtotTrkNorm->Fill(hcal_etottracknorm);
		      H_hHodBetaNtrk->Fill(hhod_beta_ntrk);
		      H_hHodBetaTrk->Fill(hhod_gtr_beta);

		      //Fill SHMS Detectors
		      H_pNGCerNpeSum->Fill(pngcer_npesum);
		      H_pHGCerNpeSum->Fill(phgcer_npesum);
		      H_pAeroNpeSum->Fill(paero_npesum);
		      H_pCalEtotNorm->Fill(pcal_etotnorm);
		      H_pCalEtotTrkNorm->Fill(pcal_etottracknorm);
		      H_pHodBetaNtrk->Fill(phod_beta_ntrk);
		      H_pHodBetaTrk->Fill(phod_gtr_beta);

		      //Fill 2D PID Correlations
		      H_hcal_vs_hcer->Fill(hcal_etottracknorm, hcer_npesum);
		      H_pcal_vs_phgcer->Fill(pcal_etottracknorm, phgcer_npesum);  
		      H_pcal_vs_pngcer->Fill(pcal_etottracknorm, pngcer_npesum);  
		      H_pcal_vs_paero->Fill(pcal_etottracknorm, paero_npesum);   
		      H_paero_vs_phgcer->Fill(paero_npesum, phgcer_npesum); 
		      H_paero_vs_pngcer->Fill(paero_npesum, pngcer_npesum); 
		      H_pngcer_vs_phgcer->Fill(pngcer_npesum, phgcer_npesum);

		      

		      //--------------------------------------------------------
		      //---------HISTOGRAM CATEGORY: Kinematics  (KIN)----------
		      //--------------------------------------------------------

		      //Fill Primary Kin Histos
		      H_the    ->Fill(th_e/dtr);
		      H_kf     ->Fill(kf);
	              H_W      ->Fill(W);
		      H_W2     ->Fill(W2);
		      H_Q2     ->Fill(Q2);
		      H_xbj    ->Fill(X);
		      H_nu     ->Fill(nu);
		      H_q      ->Fill(q);
		      H_qx     ->Fill(qx);
		      H_qy     ->Fill(qy);
		      H_qz     ->Fill(qz);
		      H_thq    ->Fill(th_q/dtr);
		      H_phq    ->Fill(ph_q/dtr);
		      H_epsilon->Fill(epsilon); 
		      
		      //Fill Secondary Kin Histos
		      H_Em       ->Fill(Em);
		      H_Em_nuc   ->Fill(Em_nuc);
		      H_Pm       ->Fill(Pm);
		      H_Pmx_lab  ->Fill(Pmx_lab);
		      H_Pmy_lab  ->Fill(Pmy_lab);
		      H_Pmz_lab  ->Fill(Pmz_lab);
		      H_Pmx_q    ->Fill(Pmx_q);
		      H_Pmy_q    ->Fill(Pmy_q);
		      H_Pmz_q    ->Fill(Pmz_q);
		      H_Tx       ->Fill(Tx);
		      H_Tr       ->Fill(Tr);
		      H_MM       ->Fill(MM);
		      H_MM2      ->Fill(MM2);
		      H_thx      ->Fill(th_x/dtr);
		      H_Pf       ->Fill(Pf);
		      H_thxq     ->Fill(th_xq/dtr);
		      H_thrq     ->Fill(th_rq/dtr);
		      H_phxq     ->Fill(ph_xq/dtr);
		      H_phrq     ->Fill(ph_rq/dtr);
		      H_Tx_cm    ->Fill(Tx_cm);
		      H_Tr_cm    ->Fill(Tr_cm);
		      H_thxq_cm  ->Fill(th_xq_cm/dtr);
		      H_thrq_cm  ->Fill(th_rq_cm/dtr);
		      H_phxq_cm  ->Fill(ph_xq_cm/dtr);
		      H_phrq_cm  ->Fill(ph_rq_cm/dtr);
		      H_Ttot_cm  ->Fill(Ttot_cm);
		      H_MandelS  ->Fill(MandelS);
		      H_MandelT  ->Fill(MandelT);
		      H_MandelU  ->Fill(MandelU);

		      //Fill (cosine, sine) of angles relative to q		      
		      H_cth_xq  ->Fill(cos(th_xq));
		      H_cth_rq  ->Fill(cos(th_rq));
		      H_sth_xq  ->Fill(sin(th_xq));
		      H_sth_rq  ->Fill(sin(th_rq));
		      H_cphi_xq ->Fill(cos(ph_xq));
		      H_cphi_rq ->Fill(cos(ph_rq));
		      H_sphi_xq ->Fill(sin(ph_xq));
		      H_sphi_rq ->Fill(sin(ph_rq));
		      //CM Frame
		      H_cth_xq_cm  ->Fill(cos(th_xq_cm));
		      H_cth_rq_cm  ->Fill(cos(th_rq_cm));
		      H_sth_xq_cm  ->Fill(sin(th_xq_cm));
		      H_sth_rq_cm  ->Fill(sin(th_rq_cm));
		      H_cphi_xq_cm ->Fill(cos(ph_xq_cm));
		      H_cphi_rq_cm ->Fill(cos(ph_rq_cm));
		      H_sphi_xq_cm ->Fill(sin(ph_xq_cm));
		      H_sphi_rq_cm ->Fill(sin(ph_rq_cm));

		  

		      //----------------------------------------------------------------------
		      //---------HISTOGRAM CATEGORY: Spectrometer Acceptance  (ACCP)----------
		      //----------------------------------------------------------------------
		      //Fill SPECTROMETER  ACCEPTANCE
		      H_exfp       ->Fill(e_xfp);
		      H_eyfp       ->Fill(e_yfp);
		      H_expfp      ->Fill(e_xpfp);
		      H_eypfp      ->Fill(e_ypfp);
		      
		      H_eytar      ->Fill(e_ytar);
		      H_exptar     ->Fill(e_xptar);
		      H_eyptar     ->Fill(e_yptar);
		      H_edelta     ->Fill(e_delta);
		      
		      H_hxfp       ->Fill(h_xfp);
		      H_hyfp       ->Fill(h_yfp);
		      H_hxpfp      ->Fill(h_xpfp);
		      H_hypfp      ->Fill(h_ypfp);
		      
		      H_hytar       ->Fill(h_ytar);
		      H_hxptar      ->Fill(h_xptar);
		      H_hyptar      ->Fill(h_yptar);
		      H_hdelta      ->Fill(h_delta);
		      
		      H_htar_x       ->Fill(htar_x);
		      H_htar_y       ->Fill(htar_y);
		      H_htar_z       ->Fill(htar_z);
		      H_etar_x       ->Fill(etar_x);
		      H_etar_y       ->Fill(etar_y);
		      H_etar_z       ->Fill(etar_z);
		      H_ztar_diff    ->Fill(ztar_diff);
		      
		      H_hXColl      ->Fill(hXColl);
		      H_hYColl      ->Fill(hYColl);
		      H_eXColl      ->Fill(eXColl);
		      H_eYColl      ->Fill(eYColl);
		      
		      H_hXColl_vs_hYColl  ->Fill(hYColl, hXColl);
		      H_eXColl_vs_eYColl  ->Fill(eYColl, eXColl);
		      
		      H_hxfp_vs_hyfp  ->Fill(h_yfp, h_xfp);
		      H_exfp_vs_eyfp  ->Fill(e_yfp, e_xfp);
		      
		      
		    }
		  //----------------------END: Fill DATA Histograms-----------------------
		  
		  
		}

	      //------END: REQUIRE "NO EDTM" CUT TO FILL DATA HISTOGRAMS-----
	      
	    }

	  //-----END: BCM Current Cut------

	  
	  //Increment Scaler Read if event == scaler_evt_perlimit for that scaler read
	  //Explanation: each scaler read value has an associated upper limit in the event number
	  //(i.e., 0->233, 1->1231, 2->2455, . . .), therefore, when we loop over data even number (gevnum)
	  // we check it against the upper limit (scal_evt_num[scal_read]) event number of a particular scaler read
	  //Once the data event number reaches that upper limit, then we increment the scaler read by +1, to go to the
	  //next scaler read, which has a new data event upper limit
	  // if there is no gevnum, then this can be replaced by (ientry + 1)
	  //cout << "gevnum = " << std::setprecision(5) << gevnum << " " << " ientry = " << (ientry + 1) << " diff = " << (ientry + 1) - gevnum << endl;
	  
	  //gevnum = ientry + 1;

	  if(gevnum==scal_evt_num[scal_read]){ scal_read++; }
	  

	  cout << "DataEventLoop: " << std::setprecision(2) << double(ientry) / nentries * 100. << "  % " << std::flush << "\r";

	}//END DATA EVENT LOOP

    }//END DATA ANALYSIS

  if(analysis=="simc")
    {
      cout << "SIMC ANALYSIS needs to be done . . . " << endl;
    }
  
}

//_______________________________________________________________________________
void baseAnalyzer::CalcEff()
{
  cout << "Calling Base CalcEff() . . . " << endl;

  //Brief: In this method, the total charge, live time, tracking efficiencies are calculated

  //Calculate Average BCM Current                                                  
  avg_current_bcm_cut = total_charge_bcm_cut / total_time_bcm_cut; //uA                              
  
  //Convert charge from uC to mC                                   
  total_charge_bcm_cut = total_charge_bcm_cut / 1000.; 

  //Convert Trigger/EDTM Rates from Hz to kHz 
  S1XscalerRate_bcm_cut   = S1XscalerRate_bcm_cut   / 1000.;
  TRIG1scalerRate_bcm_cut = TRIG1scalerRate_bcm_cut / 1000.;
  TRIG2scalerRate_bcm_cut = TRIG2scalerRate_bcm_cut / 1000.;
  TRIG3scalerRate_bcm_cut = TRIG3scalerRate_bcm_cut / 1000.;
  TRIG4scalerRate_bcm_cut = TRIG4scalerRate_bcm_cut / 1000.;
  TRIG5scalerRate_bcm_cut = TRIG5scalerRate_bcm_cut / 1000.;
  TRIG6scalerRate_bcm_cut = TRIG6scalerRate_bcm_cut / 1000.;
  EDTMscalerRate_bcm_cut  = EDTMscalerRate_bcm_cut  / 1000.;

  //Calculate Pure Computer Live Time (numerator->accepted tdc trig requires NO EDTM :: denominator -> EDTM has already been subtracted from scaler counts)
  //Pre-Scale factor has been accounted 
  cpuLT_trig1 = total_trig1_accp_bcm_cut / (total_trig1_scaler_bcm_cut / Ps1_factor);
  cpuLT_trig2 = total_trig2_accp_bcm_cut / (total_trig2_scaler_bcm_cut / Ps2_factor);
  cpuLT_trig3 = total_trig3_accp_bcm_cut / (total_trig3_scaler_bcm_cut / Ps3_factor);
  cpuLT_trig4 = total_trig4_accp_bcm_cut / (total_trig4_scaler_bcm_cut / Ps4_factor);
  cpuLT_trig5 = total_trig5_accp_bcm_cut / (total_trig5_scaler_bcm_cut / Ps5_factor);
  cpuLT_trig6 = total_trig6_accp_bcm_cut / (total_trig6_scaler_bcm_cut / Ps6_factor);

  //Calculate Computer Live Time Error (Use Binomial Statistics Error formula: err^2 = N * P * (1-P), where S->total counts (or trials), and P->probability of success : accepted_triggers  / scalers 
  cpuLT_trig1_err_Bi = sqrt( total_trig1_accp_bcm_cut * (1. - (total_trig1_accp_bcm_cut )/total_trig1_scaler_bcm_cut ) ) * Ps1_factor / total_trig1_scaler_bcm_cut;
  cpuLT_trig2_err_Bi = sqrt( total_trig2_accp_bcm_cut * (1. - (total_trig2_accp_bcm_cut )/total_trig2_scaler_bcm_cut ) ) * Ps2_factor / total_trig2_scaler_bcm_cut;
  cpuLT_trig3_err_Bi = sqrt( total_trig3_accp_bcm_cut * (1. - (total_trig3_accp_bcm_cut )/total_trig3_scaler_bcm_cut ) ) * Ps3_factor / total_trig3_scaler_bcm_cut;
  cpuLT_trig4_err_Bi = sqrt( total_trig4_accp_bcm_cut * (1. - (total_trig4_accp_bcm_cut )/total_trig4_scaler_bcm_cut ) ) * Ps4_factor / total_trig4_scaler_bcm_cut;
  cpuLT_trig5_err_Bi = sqrt( total_trig5_accp_bcm_cut * (1. - (total_trig5_accp_bcm_cut )/total_trig5_scaler_bcm_cut ) ) * Ps5_factor / total_trig5_scaler_bcm_cut;
  cpuLT_trig6_err_Bi = sqrt( total_trig6_accp_bcm_cut * (1. - (total_trig6_accp_bcm_cut )/total_trig6_scaler_bcm_cut ) ) * Ps6_factor / total_trig6_scaler_bcm_cut;

  //Calculate Computer Live Time Error (Use Bayesian Statistics Error formula)
  cpuLT_trig1_err_Bay = sqrt( (total_trig1_accp_bcm_cut * Ps1_factor + 1)*(total_trig1_accp_bcm_cut * Ps1_factor + 2)/((total_trig1_scaler_bcm_cut + 2)*(total_trig1_scaler_bcm_cut + 3)) - pow((total_trig1_accp_bcm_cut * Ps1_factor + 1),2)/pow((total_trig1_scaler_bcm_cut + 2),2) );
  cpuLT_trig2_err_Bay = sqrt( (total_trig2_accp_bcm_cut * Ps2_factor + 1)*(total_trig2_accp_bcm_cut * Ps2_factor + 2)/((total_trig2_scaler_bcm_cut + 2)*(total_trig2_scaler_bcm_cut + 3)) - pow((total_trig2_accp_bcm_cut * Ps2_factor + 1),2)/pow((total_trig2_scaler_bcm_cut + 2),2) );
  cpuLT_trig3_err_Bay = sqrt( (total_trig3_accp_bcm_cut * Ps3_factor + 1)*(total_trig3_accp_bcm_cut * Ps3_factor + 2)/((total_trig3_scaler_bcm_cut + 2)*(total_trig3_scaler_bcm_cut + 3)) - pow((total_trig3_accp_bcm_cut * Ps3_factor + 1),2)/pow((total_trig3_scaler_bcm_cut + 2),2) );
  cpuLT_trig4_err_Bay = sqrt( (total_trig4_accp_bcm_cut * Ps4_factor + 1)*(total_trig4_accp_bcm_cut * Ps4_factor + 2)/((total_trig4_scaler_bcm_cut + 2)*(total_trig4_scaler_bcm_cut + 3)) - pow((total_trig4_accp_bcm_cut * Ps4_factor + 1),2)/pow((total_trig4_scaler_bcm_cut + 2),2) );
  cpuLT_trig5_err_Bay = sqrt( (total_trig5_accp_bcm_cut * Ps5_factor + 1)*(total_trig5_accp_bcm_cut * Ps5_factor + 2)/((total_trig5_scaler_bcm_cut + 2)*(total_trig5_scaler_bcm_cut + 3)) - pow((total_trig5_accp_bcm_cut * Ps5_factor + 1),2)/pow((total_trig5_scaler_bcm_cut + 2),2) );
  cpuLT_trig6_err_Bay = sqrt( (total_trig6_accp_bcm_cut * Ps6_factor + 1)*(total_trig6_accp_bcm_cut * Ps6_factor + 2)/((total_trig6_scaler_bcm_cut + 2)*(total_trig6_scaler_bcm_cut + 3)) - pow((total_trig6_accp_bcm_cut * Ps6_factor + 1),2)/pow((total_trig6_scaler_bcm_cut + 2),2) );
  
  //Calculated total EDTM Live Time
  tLT_trig1 = total_edtm_accp_bcm_cut / (total_edtm_scaler_bcm_cut / Ps1_factor);  
  tLT_trig2 = total_edtm_accp_bcm_cut / (total_edtm_scaler_bcm_cut / Ps2_factor);  
  tLT_trig3 = total_edtm_accp_bcm_cut / (total_edtm_scaler_bcm_cut / Ps3_factor);  
  tLT_trig4 = total_edtm_accp_bcm_cut / (total_edtm_scaler_bcm_cut / Ps4_factor);  
  tLT_trig5 = total_edtm_accp_bcm_cut / (total_edtm_scaler_bcm_cut / Ps5_factor);  
  tLT_trig6 = total_edtm_accp_bcm_cut / (total_edtm_scaler_bcm_cut / Ps6_factor);  

  //Calculate EDTM Live Time Error (Use Binomial Error)
  tLT_trig1_err_Bi = sqrt( total_edtm_accp_bcm_cut * (1. - (total_edtm_accp_bcm_cut )/total_edtm_scaler_bcm_cut ) ) * Ps1_factor / total_edtm_scaler_bcm_cut;
  tLT_trig2_err_Bi = sqrt( total_edtm_accp_bcm_cut * (1. - (total_edtm_accp_bcm_cut )/total_edtm_scaler_bcm_cut ) ) * Ps2_factor / total_edtm_scaler_bcm_cut;
  tLT_trig3_err_Bi = sqrt( total_edtm_accp_bcm_cut * (1. - (total_edtm_accp_bcm_cut )/total_edtm_scaler_bcm_cut ) ) * Ps3_factor / total_edtm_scaler_bcm_cut;
  tLT_trig4_err_Bi = sqrt( total_edtm_accp_bcm_cut * (1. - (total_edtm_accp_bcm_cut )/total_edtm_scaler_bcm_cut ) ) * Ps4_factor / total_edtm_scaler_bcm_cut;
  tLT_trig5_err_Bi = sqrt( total_edtm_accp_bcm_cut * (1. - (total_edtm_accp_bcm_cut )/total_edtm_scaler_bcm_cut ) ) * Ps5_factor / total_edtm_scaler_bcm_cut;
  tLT_trig6_err_Bi = sqrt( total_edtm_accp_bcm_cut * (1. - (total_edtm_accp_bcm_cut )/total_edtm_scaler_bcm_cut ) ) * Ps6_factor / total_edtm_scaler_bcm_cut;

  //Calculate EDTM Live Time Error (Use Bayesian Error)
  tLT_trig1_err_Bay = sqrt( (total_edtm_accp_bcm_cut * Ps1_factor + 1)*(total_edtm_accp_bcm_cut * Ps1_factor + 2)/((total_edtm_scaler_bcm_cut + 2)*(total_edtm_scaler_bcm_cut + 3)) - pow((total_edtm_accp_bcm_cut * Ps1_factor + 1),2)/pow((total_edtm_scaler_bcm_cut + 2),2) );
  tLT_trig2_err_Bay = sqrt( (total_edtm_accp_bcm_cut * Ps2_factor + 1)*(total_edtm_accp_bcm_cut * Ps2_factor + 2)/((total_edtm_scaler_bcm_cut + 2)*(total_edtm_scaler_bcm_cut + 3)) - pow((total_edtm_accp_bcm_cut * Ps2_factor + 1),2)/pow((total_edtm_scaler_bcm_cut + 2),2) );
  tLT_trig3_err_Bay = sqrt( (total_edtm_accp_bcm_cut * Ps3_factor + 1)*(total_edtm_accp_bcm_cut * Ps3_factor + 2)/((total_edtm_scaler_bcm_cut + 2)*(total_edtm_scaler_bcm_cut + 3)) - pow((total_edtm_accp_bcm_cut * Ps3_factor + 1),2)/pow((total_edtm_scaler_bcm_cut + 2),2) );
  tLT_trig4_err_Bay = sqrt( (total_edtm_accp_bcm_cut * Ps4_factor + 1)*(total_edtm_accp_bcm_cut * Ps4_factor + 2)/((total_edtm_scaler_bcm_cut + 2)*(total_edtm_scaler_bcm_cut + 3)) - pow((total_edtm_accp_bcm_cut * Ps4_factor + 1),2)/pow((total_edtm_scaler_bcm_cut + 2),2) );
  tLT_trig5_err_Bay = sqrt( (total_edtm_accp_bcm_cut * Ps5_factor + 1)*(total_edtm_accp_bcm_cut * Ps5_factor + 2)/((total_edtm_scaler_bcm_cut + 2)*(total_edtm_scaler_bcm_cut + 3)) - pow((total_edtm_accp_bcm_cut * Ps5_factor + 1),2)/pow((total_edtm_scaler_bcm_cut + 2),2) );
  tLT_trig6_err_Bay = sqrt( (total_edtm_accp_bcm_cut * Ps6_factor + 1)*(total_edtm_accp_bcm_cut * Ps6_factor + 2)/((total_edtm_scaler_bcm_cut + 2)*(total_edtm_scaler_bcm_cut + 3)) - pow((total_edtm_accp_bcm_cut * Ps6_factor + 1),2)/pow((total_edtm_scaler_bcm_cut + 2),2) );
  
  //Ensure that if Ps_factor = -1 (trigger input OFF), then live times default to -1.0
  if(Ps1_factor==-1) { cpuLT_trig1 = -1.0, tLT_trig1 = -1.0; }
  if(Ps2_factor==-1) { cpuLT_trig2 = -1.0, tLT_trig2 = -1.0; }
  if(Ps3_factor==-1) { cpuLT_trig3 = -1.0, tLT_trig3 = -1.0; }
  if(Ps4_factor==-1) { cpuLT_trig4 = -1.0, tLT_trig4 = -1.0; }
  if(Ps5_factor==-1) { cpuLT_trig5 = -1.0, tLT_trig5 = -1.0; }
  if(Ps6_factor==-1) { cpuLT_trig6 = -1.0, tLT_trig6 = -1.0; }
  
  //Choose what trigger type to use in correction factor (see set_basic_cuts.inp)
  if(trig_type=="trig1") { trig_rate = TRIG1scalerRate_bcm_cut, cpuLT_trig = cpuLT_trig1, cpuLT_trig_err_Bi = cpuLT_trig1_err_Bi, cpuLT_trig_err_Bay = cpuLT_trig1_err_Bay, tLT_trig = tLT_trig1, tLT_trig_err_Bi = tLT_trig1_err_Bi, tLT_trig_err_Bay = tLT_trig1_err_Bay, Ps_factor = Ps1_factor, total_trig_scaler_bcm_cut = total_trig1_scaler_bcm_cut, total_trig_accp_bcm_cut = total_trig1_accp_bcm_cut; }
  if(trig_type=="trig2") { trig_rate = TRIG2scalerRate_bcm_cut, cpuLT_trig = cpuLT_trig2, cpuLT_trig_err_Bi = cpuLT_trig2_err_Bi, cpuLT_trig_err_Bay = cpuLT_trig2_err_Bay, tLT_trig = tLT_trig2, tLT_trig_err_Bi = tLT_trig2_err_Bi, tLT_trig_err_Bay = tLT_trig2_err_Bay, Ps_factor = Ps2_factor, total_trig_scaler_bcm_cut = total_trig2_scaler_bcm_cut, total_trig_accp_bcm_cut = total_trig2_accp_bcm_cut; }
  if(trig_type=="trig3") { trig_rate = TRIG3scalerRate_bcm_cut, cpuLT_trig = cpuLT_trig3, cpuLT_trig_err_Bi = cpuLT_trig3_err_Bi, cpuLT_trig_err_Bay = cpuLT_trig3_err_Bay, tLT_trig = tLT_trig3, tLT_trig_err_Bi = tLT_trig3_err_Bi, tLT_trig_err_Bay = tLT_trig3_err_Bay, Ps_factor = Ps3_factor, total_trig_scaler_bcm_cut = total_trig3_scaler_bcm_cut, total_trig_accp_bcm_cut = total_trig3_accp_bcm_cut; }
  if(trig_type=="trig4") { trig_rate = TRIG4scalerRate_bcm_cut, cpuLT_trig = cpuLT_trig4, cpuLT_trig_err_Bi = cpuLT_trig4_err_Bi, cpuLT_trig_err_Bay = cpuLT_trig4_err_Bay, tLT_trig = tLT_trig4, tLT_trig_err_Bi = tLT_trig4_err_Bi, tLT_trig_err_Bay = tLT_trig4_err_Bay, Ps_factor = Ps4_factor, total_trig_scaler_bcm_cut = total_trig4_scaler_bcm_cut, total_trig_accp_bcm_cut = total_trig4_accp_bcm_cut; }
  if(trig_type=="trig5") { trig_rate = TRIG5scalerRate_bcm_cut, cpuLT_trig = cpuLT_trig5, cpuLT_trig_err_Bi = cpuLT_trig5_err_Bi, cpuLT_trig_err_Bay = cpuLT_trig5_err_Bay, tLT_trig = tLT_trig5, tLT_trig_err_Bi = tLT_trig5_err_Bi, tLT_trig_err_Bay = tLT_trig5_err_Bay, Ps_factor = Ps5_factor, total_trig_scaler_bcm_cut = total_trig5_scaler_bcm_cut, total_trig_accp_bcm_cut = total_trig5_accp_bcm_cut; }
  if(trig_type=="trig6") { trig_rate = TRIG6scalerRate_bcm_cut, cpuLT_trig = cpuLT_trig6, cpuLT_trig_err_Bi = cpuLT_trig6_err_Bi, cpuLT_trig_err_Bay = cpuLT_trig6_err_Bay, tLT_trig = tLT_trig6, tLT_trig_err_Bi = tLT_trig6_err_Bi, tLT_trig_err_Bay = tLT_trig6_err_Bay, Ps_factor = Ps6_factor, total_trig_scaler_bcm_cut = total_trig6_scaler_bcm_cut, total_trig_accp_bcm_cut = total_trig6_accp_bcm_cut; }


  //Calculate HMS Tracking Efficiency                                                                                                                 
  hTrkEff = h_did / h_should;                                                                                                                  
  hTrkEff_err = sqrt(h_should-h_did) / h_should;
  
  //Calculate SHMS Tracking Efficiency                                                                                                
  pTrkEff = p_did / p_should; 
  pTrkEff_err = sqrt(p_should-p_did) / p_should;                                                            
  
    
}

//_______________________________________________________________________________
void baseAnalyzer::ApplyWeight()
{
  /*
    Brief: Method to apply efficiency corrections to the data Yield on a run-by-run basis, 
    by scaling the Filled Histograms. Corrections applied include : 
   
    1) normalization by total_charge 
    2) total (EDTM) live time correction, 
    3) tracking efficiency correction (HMS, SHMS), 
    4) hadron (proton, pion or Kaon) absorption correction, 
    5) target boiling corrections 
    
    NOTE: Additional corrections may be added, depending on the experimental analysis
  */
  cout << "Calling Base ApplyWeight() . . . " << endl; 

  /*
    Target Boiling Slopes (May change depending on the target used)
    Currently, these slopes are from boiling studies done on April 02, 2018,
    However, these must be updated with boiling studies results done during other experiments
    
    //Target Boiling Slopes/Error (C. Yero Boiling Studies from April 02, 2018) 
    //Units (fractional Yield_loss / uA)     
    
    LH2_slope = 6.34e-4,  LH2_slope_err = 6.2e-5;  
    LD2_slope = 8.00e-4;  LD2_slope_err = 7.0e-5;

    dIb_Ib = 0.02;  //For now, assume that the BCM relative error is 2% (0.02)
    
    //Absolute error on the average beam current is:
    avg_current_bcm_cut_err = dIb_Ib * avg_current_bcm_cut;

    //example: this correction represents fraction of events NOT lost due to boiling, 
    //e.g. 0.98 -> 98% of events NOT lost due to boiling. Dividing the Yield by 0.98 to recover the 2% of the events lost (yield count increases)
    tgtBoil_corr = (1. - LH2_slope * avg_current_bcm_cut)
    tgtBoil_corr_err = sqrt( pow(avg_current_bcm_cut*LH2_slope_err ,2) + pow(LH2_slope*avg_current_bcm_cut_err ,2) )
    
  */
  
  //For now, assume that the BCM relative error is 2% (0.02)
  dIb_Ib = 0.02;  

  //Absolute error on the average beam current is:
  avg_current_bcm_cut_err = dIb_Ib * avg_current_bcm_cut;
  
  //April 02, 2018 (C. Yero tgt. boiling study results for HMS)
  LH2_slope = 6.34e-4,  LH2_slope_err = 6.2e-5;  
  LD2_slope = 8.00e-4;  LD2_slope_err = 7.0e-5;
  

  if(tgt_type=="LH2") { tgt_slope = LH2_slope, tgt_slope_err = LH2_slope_err;}
  if(tgt_type=="LD2") { tgt_slope = LD2_slope, tgt_slope_err = LD2_slope_err;}

  tgtBoil_corr = (1. - tgt_slope * avg_current_bcm_cut);
  tgtBoil_corr_err = sqrt( pow( avg_current_bcm_cut * tgt_slope_err , 2)  + pow( tgt_slope * avg_current_bcm_cut_err ,2) );

  //For now, assume no target boiling
  tgtBoil_corr = 1.;
  tgtBoil_corr_err = -1.0;
  
  /*
    Proton, Pion or Kaon absorption correction:

    This correction factor is calculated by estimating the fraction of
    hadrons that actually made it to the hodoscopoes to form a trigger 
    divided by the number of hadrons that should have made it to form a trigger.
    This calculation can be thought of an efficiency calculation for hadrons that
    were not absorbed by the materials (target, entrance/exit windows, detector material, etc.)

    Depending on the nature of the experiment, other results for other hadron absorption 
    corrections can be used.
  
    //hadAbs_corr = 1.;         (For now assume there is no hadrons absorbed) 
    //hadAbs_corr_err = -1.;    (For now assume there is no hadrons absorbed)
    
    //C. Yero proton absorption results for HMS during E12-20-003 Commissioning
    hadAbs_corr     = 0.9534;
    hadAbs_corr_err = 0.0047; 

  */
  
  //C. Yero proton absorption results for HMS during E12-20-003 Commissioning
  hadAbs_corr     = 0.9534;
  hadAbs_corr_err = 0.0047; 

  //For now, assume no hadron absorption
  hadAbs_corr     = 1.;
  hadAbs_corr_err = 1.; 
  
  
  //Full Weight
  //FullWeight = 1. / (total_charge_bcm_cut * hTrkEff * pTrkEff * tLT_trig * tgtBoil_corr * hadAbs_corr);

  //For testing purposes, only normalize by total charge
  FullWeight = 1.; // / total_charge_bcm_cut;
  
  //Scale Data Histograms by Full Weight (Each run for a particular kinematics can then be combined, once they are scaled by the FullWeight)

  //--------------------------------------------------------------------
  //---------HISTOGRAM CATEGORY: Particle Identification (PID)----------
  //--------------------------------------------------------------------

  //Scale Coincidence Time		      
  H_ep_ctime->Scale(FullWeight);
  H_eK_ctime->Scale(FullWeight);
  H_ePi_ctime->Scale(FullWeight);
  
  //Scale HMS Detectors
  H_hCerNpeSum->Scale(FullWeight);
  H_hCalEtotNorm->Scale(FullWeight);
  H_hCalEtotTrkNorm->Scale(FullWeight);
  H_hHodBetaNtrk->Scale(FullWeight);
  H_hHodBetaTrk->Scale(FullWeight);
  
  //Scale SHMS Detectors
  H_pNGCerNpeSum->Scale(FullWeight);
  H_pHGCerNpeSum->Scale(FullWeight);
  H_pAeroNpeSum->Scale(FullWeight);
  H_pCalEtotNorm->Scale(FullWeight);
  H_pCalEtotTrkNorm->Scale(FullWeight);
  H_pHodBetaNtrk->Scale(FullWeight);
  H_pHodBetaTrk->Scale(FullWeight);

  //Scale 2D  PID Correlations
  H_hcal_vs_hcer->Scale(FullWeight);
  H_pcal_vs_phgcer->Scale(FullWeight);  
  H_pcal_vs_pngcer->Scale(FullWeight);  
  H_pcal_vs_paero->Scale(FullWeight);   
  H_paero_vs_phgcer->Scale(FullWeight); 
  H_paero_vs_pngcer->Scale(FullWeight); 
  H_pngcer_vs_phgcer->Scale(FullWeight);
  
  //--------------------------------------------------------
  //---------HISTOGRAM CATEGORY: Kinematics  (KIN)----------
  //--------------------------------------------------------
  
  //Fill Primary Kin Histos
  H_the    ->Scale(FullWeight);
  H_kf     ->Scale(FullWeight);
  H_W      ->Scale(FullWeight);
  H_W2     ->Scale(FullWeight);
  H_Q2     ->Scale(FullWeight);
  H_xbj    ->Scale(FullWeight);
  H_nu     ->Scale(FullWeight);
  H_q      ->Scale(FullWeight);
  H_qx     ->Scale(FullWeight);
  H_qy     ->Scale(FullWeight);
  H_qz     ->Scale(FullWeight);
  H_thq    ->Scale(FullWeight);
  H_phq    ->Scale(FullWeight);
  H_epsilon->Scale(FullWeight);
  
  //Fill Secondary Kin Histos
  H_Em       ->Scale(FullWeight);
  H_Em_nuc   ->Scale(FullWeight);
  H_Pm       ->Scale(FullWeight);
  H_Pmx_lab  ->Scale(FullWeight);
  H_Pmy_lab  ->Scale(FullWeight);
  H_Pmz_lab  ->Scale(FullWeight);
  H_Pmx_q    ->Scale(FullWeight);
  H_Pmy_q    ->Scale(FullWeight);
  H_Pmz_q    ->Scale(FullWeight);
  H_Tx       ->Scale(FullWeight);
  H_Tr       ->Scale(FullWeight);
  H_MM       ->Scale(FullWeight);
  H_MM2      ->Scale(FullWeight);
  H_thx      ->Scale(FullWeight);
  H_Pf       ->Scale(FullWeight);
  H_thxq     ->Scale(FullWeight);
  H_thrq     ->Scale(FullWeight);
  H_phxq     ->Scale(FullWeight);
  H_phrq     ->Scale(FullWeight);
  H_Tx_cm    ->Scale(FullWeight);
  H_Tr_cm    ->Scale(FullWeight);
  H_thxq_cm  ->Scale(FullWeight);
  H_thrq_cm  ->Scale(FullWeight);
  H_phxq_cm  ->Scale(FullWeight);
  H_phrq_cm  ->Scale(FullWeight);
  H_Ttot_cm  ->Scale(FullWeight);
  H_MandelS  ->Scale(FullWeight);
  H_MandelT  ->Scale(FullWeight);
  H_MandelU  ->Scale(FullWeight);

  //Fill (cosine, sine) of angles relative to q		      
  H_cth_xq  ->Scale(FullWeight);
  H_cth_rq  ->Scale(FullWeight);
  H_sth_xq  ->Scale(FullWeight);
  H_sth_rq  ->Scale(FullWeight);
  H_cphi_xq ->Scale(FullWeight);
  H_cphi_rq ->Scale(FullWeight);
  H_sphi_xq ->Scale(FullWeight);
  H_sphi_rq ->Scale(FullWeight);
  //CM Frame
  H_cth_xq_cm  ->Scale(FullWeight);
  H_cth_rq_cm  ->Scale(FullWeight);
  H_sth_xq_cm  ->Scale(FullWeight);
  H_sth_rq_cm  ->Scale(FullWeight);
  H_cphi_xq_cm ->Scale(FullWeight);
  H_cphi_rq_cm ->Scale(FullWeight);
  H_sphi_xq_cm ->Scale(FullWeight);
  H_sphi_rq_cm ->Scale(FullWeight);

  

  //----------------------------------------------------------------------
  //---------HISTOGRAM CATEGORY: Spectrometer Acceptance  (ACCP)----------
  //----------------------------------------------------------------------
  
  //Add ACCP Histos to TList
  H_exfp       ->Scale(FullWeight);
  H_eyfp       ->Scale(FullWeight);
  H_expfp      ->Scale(FullWeight);
  H_eypfp      ->Scale(FullWeight);
  
  H_eytar      ->Scale(FullWeight);
  H_exptar     ->Scale(FullWeight);
  H_eyptar     ->Scale(FullWeight);
  H_edelta     ->Scale(FullWeight);
  
  H_hxfp       ->Scale(FullWeight);
  H_hyfp       ->Scale(FullWeight);
  H_hxpfp      ->Scale(FullWeight);
  H_hypfp      ->Scale(FullWeight);
  
  H_hytar       ->Scale(FullWeight);
  H_hxptar      ->Scale(FullWeight);
  H_hyptar      ->Scale(FullWeight);
  H_hdelta      ->Scale(FullWeight);
  
  H_htar_x       ->Scale(FullWeight);
  H_htar_y       ->Scale(FullWeight);
  H_htar_z       ->Scale(FullWeight);
  H_etar_x       ->Scale(FullWeight);
  H_etar_y       ->Scale(FullWeight);
  H_etar_z       ->Scale(FullWeight);
  H_ztar_diff    ->Scale(FullWeight);
  
  H_hXColl      ->Scale(FullWeight);
  H_hYColl      ->Scale(FullWeight);
  H_eXColl      ->Scale(FullWeight);
  H_eYColl      ->Scale(FullWeight);
  
  H_hXColl_vs_hYColl  ->Scale(FullWeight);
  H_eXColl_vs_eYColl  ->Scale(FullWeight);
  
  H_hxfp_vs_hyfp  ->Scale(FullWeight);
  H_exfp_vs_eyfp  ->Scale(FullWeight);
  
}

//_______________________________________________________________________________
void baseAnalyzer::WriteHist()
{
  /*
    Brief: Method to write histograms to a ROOTfile
  */

  cout << "Calling WriteHist() . . ." << endl;

  
  //Write Data Histograms
  if(analysis=="data")
    {
      //Create Output ROOTfile
      outROOT = new TFile(data_OutputFileName, "RECREATE");

      //Make directories to store histograms based on category
      outROOT->mkdir("pid_plots");
      outROOT->mkdir("kin_plots");
      outROOT->mkdir("accp_plots");

      //Write PID histos to pid_plots directory
      outROOT->cd("pid_plots");
      pid_HList->Write();
      
      //Write Kinematics histos to kin_plots directory
      outROOT->cd("kin_plots");
      kin_HList->Write();

      //Write Acceptance histos to accp_plots directory
      outROOT->cd("accp_plots");
      accp_HList->Write();

      //Close File
      outROOT->Close();
    }

  
}

//_______________________________________________________________________________
void baseAnalyzer::WriteReport()
{
  
  /*Method to write charge, efficiencies, live time and other relevant quantities to a data file*/
  
  cout << "Calling WriteReport() . . ." << endl;

  
  if(analysis=="data"){

    //---------------------------------------------------------

    //Choose which trigger is relevatn
    
    //Check if file already exists
    in_file.open(report_OutputFileName.Data());
    
    if(in_file.fail()){
      
      cout << "Report File does NOT exist, will create one . . . " << endl;
      
      out_file.open(report_OutputFileName);
      out_file << "#-------------------------------------" << endl;
      out_file << "#        Data Analysis Summary        " << endl;
      out_file << "#-------------------------------------" << endl;
      out_file << "#                                     " << endl;
      out_file << Form("# %s  | Beam Current Threshold: > %.2f uA ", bcm_type.Data(), bcm_thrs) << endl;
      out_file << "#                                     " << endl;
      out_file << Form("# DAQ Mode: %s | Trigger: %s              ", daq_mode.Data(), trig_type.Data()) << endl;
      out_file << Form("# electron arm: %s                        ", e_arm_name.Data() ) << endl;
      out_file << "#                                              " << endl;
      out_file << "#---PID Cuts--- " << endl;
      if(eKctime_pidCut_flag)         {out_file << Form("# Kaon Coincidence Time Cut: (%.3f,%.3f) ns", cpid_eKctime_min, cpid_eKctime_max) << endl;}
      if(ePictime_pidCut_flag)        {out_file << Form("# Pion Coincidence Time Cut: (%.3f,%.3f) ns", cpid_ePictime_min, cpid_ePictime_max) << endl;}
      if(ePctime_pidCut_flag)         {out_file << Form("# Proton Coincidence Time Cut: (%.3f,%.3f) ns", cpid_ePctime_min, cpid_ePctime_max) << endl;}
      if(petot_trkNorm_pidCut_flag)   {out_file << Form("# SHMS Calorimeter EtotTrackNorm Cut: (%.3f, %.3f)", cpid_petot_trkNorm_min,  cpid_petot_trkNorm_max) << endl;}
      if(pngcer_pidCut_flag) {out_file << Form("# SHMS Noble Gas Cherenkov NPE Sum Cut: (%.3f, %.3f)", cpid_pngcer_npeSum_min,  cpid_pngcer_npeSum_max) << endl;}
      if(phgcer_pidCut_flag) {out_file << Form("# SHMS Heavy Gas Cherenkov NPE Sum Cut: (%.3f, %.3f)", cpid_phgcer_npeSum_min,  cpid_phgcer_npeSum_max) << endl;}
      if(paero_pidCut_flag) {out_file << Form("# SHMS Aerogel Cherenkov NPE Sum Cut: (%.3f, %.3f)", cpid_paero_npeSum_min,  cpid_paero_npeSum_max) << endl;}
      if(hetot_trkNorm_pidCut_flag) {out_file << Form("# HMS Calorimeter EtotTrackNorm Cut: (%.3f, %.3f)", cpid_hetot_trkNorm_min,  cpid_hetot_trkNorm_max) << endl;}
      if(hcer_pidCut_flag) {out_file << Form("# HMS Gas Cherenkov NPE Sum Cut: (%.3f, %.3f)", cpid_hcer_npeSum_min,  cpid_hcer_npeSum_max) << endl;}      
      out_file << "#                                     " << endl;
      out_file << "#---Kinematics Cuts--- " << endl;
      if(Q2_cut_flag)            {out_file << Form("# 4-Momentum Transfer (Q^2): (%.3f, %.3f) GeV^2", c_Q2_min, c_Q2_max ) << endl;}
      if(Em_cut_flag)            {out_file << Form("# Missing Energy, Em: (%.3f, %.3f) GeV",   c_Em_min, c_Em_max ) << endl;}
      if(W_cut_flag)             {out_file << Form("# Invariant Mass, W: (%.3f, %.3f) GeV",   c_W_min,  c_W_max  ) << endl;}
      if(MM_K_cut_flag)          {out_file << Form("# Kaon Missing Mass, MM_K: (%.3f, %.3f) GeV",   c_MM_K_min,  c_MM_K_max  ) << endl;}
      if(MM_Pi_cut_flag)         {out_file << Form("# Pion Missing Mass, MM_Pi: (%.3f, %.3f) GeV",   c_MM_Pi_min,  c_MM_Pi_max  ) << endl;}
      if(MM_P_cut_flag)          {out_file << Form("# Proton Missing Mass, MM_P: (%.3f, %.3f) GeV",   c_MM_P_min,  c_MM_P_max  ) << endl;}
      out_file << "#                                     " << endl;
      out_file << "#---Acceptance Cuts--- " << endl;
      if(hdelta_cut_flag)        {out_file << Form("# HMS Momentum Acceptance: (%.3f, %.3f) %%",  c_hdelta_min, c_hdelta_max ) << endl;}
      if(edelta_cut_flag)        {out_file << Form("# SHMS Momentum Acceptance: (%.3f, %.3f) %%", c_edelta_min, c_edelta_max ) << endl;}
      if(ztarDiff_cut_flag)      {out_file << Form("# Z-Reaction Vertex Difference: (%.3f, %.3f) %%", c_ztarDiff_min, c_ztarDiff_max ) << endl;}
      out_file << "#                       " << endl;
      out_file << "# Units: charge [mC] | currnet [uA] | rates [kHz] |  efficiencies [fractional form]                       " << endl;
      out_file << "#                       " << endl;
      out_file << std::setw(2) << "#! Run[i,0]/" << std::setw(25) << "charge[f,1]/" << std::setw(25) << "avg_current[f,2]/" << std::setw(25)  << "hTrkEff[f,3]/" << std::setw(25) << "hTrkEff_err[f,4]/" << std::setw(25) << "pTrkEff[f,5]/" << std::setw(25) << "pTrkEff_err[f,6]/" << std::setw(25) << "tgt_boil_factor[f,7]/" << std::setw(30) << "tgt_boil_factor_err[f,8]" << std::setw(25) << "hadAbs_factor[f,9]/" << std::setw(30) << "hadAbs_factor_err[f,10]/" << std::setw(25) <<  "cpuLT[f,11]/" << std::setw(25) << "cpuLT_err_Bi[f,12]/" << std::setw(25) << "cpuLT_err_Bay[f,13]/" << std::setw(25) << "tLT[f,14]/" << std::setw(25) << "tLT_err_Bi[f,15]/" << std::setw(25) << "tLT_err_Bay[f,16]/" << std::setw(25) << "S1X_rate[f,17]/" << std::setw(25) << "trig_rate[f,18]/"  << std::setw(25) << "edtm_rate[f,19]/"  << std::setw(25) << "Pre_Scale[f,20]/" << std::setw(25) << "edtm_accp[f,21]/" << std::setw(25) << "edtm_scaler[f,22]/" << std::setw(25) << "trig_accp[f,23]/" << std::setw(25) << "trig_scaler[f,24]/" << endl;
            
      out_file.close();
      in_file.close();
    
    }

    //Open Report FIle in append mode
    out_file.open(report_OutputFileName, ios::out | ios::app);
    out_file << std::setw(7) << run  << std::setw(25) << total_charge_bcm_cut << std::setw(25) << avg_current_bcm_cut << std::setw(25) << hTrkEff << std::setw(25) << hTrkEff_err << std::setw(25) << pTrkEff << std::setw(25) << pTrkEff_err << std::setw(25) << tgtBoil_corr << std::setw(25) << tgtBoil_corr_err << std::setw(25) << hadAbs_corr << std::setw(25) << hadAbs_corr_err << std::setw(25) << cpuLT_trig << std::setw(25) << cpuLT_trig_err_Bi << std::setw(25) << cpuLT_trig_err_Bay << std::setw(25) << tLT_trig << std::setw(25) << tLT_trig_err_Bi << std::setw(25) << tLT_trig_err_Bay << std::setw(25) << S1XscalerRate_bcm_cut << std::setw(25) << trig_rate << std::setw(25) << EDTMscalerRate_bcm_cut << std::setw(25) << Ps_factor << std::setw(25) << total_edtm_accp_bcm_cut << std::setw(25) << (total_edtm_scaler_bcm_cut / Ps_factor) << std::setw(25) << total_trig_accp_bcm_cut << std::setw(25) << (total_trig_scaler_bcm_cut / Ps_factor) << endl;
    out_file.close();
  }
  
  cout << "Ending WriteReport() . . ." << endl;
  
} //End WriteReport()

//_______________________________________________________________________________
void baseAnalyzer::CombineHistos()
{
  /* Brief: Method to add histograms multiple runs of the same kinematics setting.
     Get the Histo for each run, open a new ROOTfile (_name_combined.root) in UPDATE mode,
     and Write the histograms to it.

     Currently, a loop over all histogram objects for each directory (pid_plots, kin_plots, accp_plots)
     is done during the summing of histograms. If a new directory with additional plots is added, make sure
     to add it to this method so that histograms may be summed.

     PENDING ISSUES: Still need to figure out more efficient way to 
     loop over each key and get number of elements. The problem is how to
     retrieve TList object from ROOTfiles and call TList obj->GetEntries().
  */


  //Decide whether to combine all histograms or NOT. (Look in main_controls.inp)
  //If the list of runs correspond to different kinematics, then they should NOT be combined (combine_runs_flag=0)
  if(combine_runs_flag==0) return;   //exit this function if combine_runs_flag == 0; 
  
  //Check if combined ROOTfile exits, otherwise, create it and add the histograms from the 1st run on the list
  Bool_t file_not_exist = gSystem->AccessPathName(data_OutputFileName_combined.Data());

  if(file_not_exist){
    cout << "Combined ROOTfile does NOT exist!!! | Will create it." << endl;

    //----making a copy of the root file does not work because the hel_directory is also copied, and so when this method gets
    //called in the derived class, it will add the helicity list a second time, for the first run, so it will double count
    //one needs to explicitly put only the lists of the base class here.

    //Copy root file from 1st run, to initialize final root file, and then keep adding to this file.
    //gSystem->CopyFile(data_OutputFileName, data_OutputFileName_combined);

    
    //Create Output ROOTfile
    outROOT = new TFile(data_OutputFileName_combined, "RECREATE");
    
    //Make directories to store histograms based on category
    outROOT->mkdir("pid_plots");
    outROOT->mkdir("kin_plots");
    outROOT->mkdir("accp_plots");
    
    //Write PID histos to pid_plots directory
    outROOT->cd("pid_plots");
    pid_HList->Write();
    
    //Write Kinematics histos to kin_plots directory
    outROOT->cd("kin_plots");
    kin_HList->Write();
    
    //Write Acceptance histos to accp_plots directory
    outROOT->cd("accp_plots");
    accp_HList->Write();

    outROOT->Close();
    
  }
  
  //If combined ROOTfile already exists, keep adding histo counts for each run
  
  else{
    
    //determine what class types are in the list
    TString class_name;
    
    //Set up histogram names/locations to find
    TString hist_name;
    TString hist_dir;
    
    
    outROOT = new TFile(data_OutputFileName_combined.Data(), "READ");  
    
    //------------------------------------------
    //Lopp over pid_HList of histogram objects
    //------------------------------------------
    for(int i=0; i<pid_HList->GetEntries(); i++)
      {
	
	//Determine object data type (as of now, either TH1F or TH2F are possible)
	class_name = pid_HList->At(i)->ClassName();
	
	//Read ith histograms in the list from current run
	if(class_name=="TH1F") {
	  //Get histogram from current run
	  h_i = (TH1F *)pid_HList->At(i);
	  hist_name = h_i->GetName();
	  //full path to histogram must be given, otherwise, histogram will NOT be retreived from ROOTfile
	  hist_dir = "pid_plots/" + hist_name;  
	  //Get cumulative histogram object stored in ROOTfile and add the h_i (histogram from ith run) to h_total (cumulative histogram) 
	  outROOT->GetObject(hist_dir, h_total); h_total->Add(h_i); outROOT->ReOpen("UPDATE"); outROOT->cd("pid_plots"); h_total->Write("", TObject::kOverwrite); outROOT->ReOpen("READ");	       	   
	  
	}
	
	if(class_name=="TH2F") {
	  //Get histogram from current run
	  h2_i = (TH2F *)pid_HList->At(i);  
	  hist_name = h2_i->GetName();
	  //full path to histogram must be given, otherwise, histogram will NOT be retreived from ROOTfile
	  hist_dir = "pid_plots/" + hist_name;  
	  //Get cumulative histogram object stored in ROOTfile and add the h_i (histogram from ith run) to h_total (cumulative histogram) 
	  outROOT->GetObject(hist_dir, h2_total); h2_total->Add(h2_i); outROOT->ReOpen("UPDATE"); outROOT->cd("pid_plots"); h2_total->Write("", TObject::kOverwrite); outROOT->ReOpen("READ");	       
	  
	}
	
      }//end loop over pid_HList

    //--------------------------------------------------------------------------------------------

    //------------------------------------------
    //Lopp over kin_HList of histogram objects
    //------------------------------------------
    for(int i=0; i<kin_HList->GetEntries(); i++)
      {
	
	//Determine object data type (as of now, either TH1F or TH2F are possible)
	class_name = kin_HList->At(i)->ClassName();
	
	//Read ith histograms in the list from current run
	if(class_name=="TH1F") {
	  //Get histogram from current run
	  h_i = (TH1F *)kin_HList->At(i);  
	  hist_name = h_i->GetName();
	  //full path to histogram must be given, otherwise, histogram will NOT be retreived from ROOTfile
	  hist_dir = "kin_plots/" + hist_name;  
	  //Get cumulative histogram object stored in ROOTfile and add the h_i (histogram from ith run) to h_total (cumulative histogram) 
	  outROOT->GetObject(hist_dir, h_total); h_total->Add(h_i); outROOT->ReOpen("UPDATE"); outROOT->cd("kin_plots"); h_total->Write("", TObject::kOverwrite); outROOT->ReOpen("READ");	       	   
	  
	}
	
	if(class_name=="TH2F") {
	  //Get histogram from current run
	  h2_i = (TH2F *)kin_HList->At(i);  
	  hist_name = h2_i->GetName();
	  //full path to histogram must be given, otherwise, histogram will NOT be retreived from ROOTfile
	  hist_dir = "kin_plots/" + hist_name;  
	  //Get cumulative histogram object stored in ROOTfile and add the h_i (histogram from ith run) to h_total (cumulative histogram) 
	  outROOT->GetObject(hist_dir, h2_total); h2_total->Add(h2_i); outROOT->ReOpen("UPDATE"); outROOT->cd("kin_plots"); h2_total->Write("", TObject::kOverwrite); outROOT->ReOpen("READ");	       
	  
	}
	
      }//end loop over kin_HList

    //--------------------------------------------------------------------------------------------

    //------------------------------------------
    //Lopp over accp_HList of histogram objects
    //------------------------------------------
    for(int i=0; i<accp_HList->GetEntries(); i++)
      {
	
	//Determine object data type (as of now, either TH1F or TH2F are possible)
	class_name = accp_HList->At(i)->ClassName();
	
	//Read ith histograms in the list from current run
	if(class_name=="TH1F") {
	  //Get histogram from current run
	  h_i = (TH1F *)accp_HList->At(i);  
	  hist_name = h_i->GetName();
	  //full path to histogram must be given, otherwise, histogram will NOT be retreived from ROOTfile
	  hist_dir = "accp_plots/" + hist_name;  
	  //Get cumulative histogram object stored in ROOTfile and add the h_i (histogram from ith run) to h_total (cumulative histogram) 
	  outROOT->GetObject(hist_dir, h_total); h_total->Add(h_i); outROOT->ReOpen("UPDATE"); outROOT->cd("accp_plots"); h_total->Write("", TObject::kOverwrite); outROOT->ReOpen("READ");	       	   
	  
	}
	
	if(class_name=="TH2F") {
	  //Get histogram from current run
	  h2_i = (TH2F *)accp_HList->At(i);
	  hist_name = h2_i->GetName();
	  //full path to histogram must be given, otherwise, histogram will NOT be retreived from ROOTfile
	  hist_dir = "accp_plots/" + hist_name;
	  //Get cumulative histogram object stored in ROOTfile and add the h_i (histogram from ith run) to h_total (cumulative histogram) 
	  outROOT->GetObject(hist_dir, h2_total); h2_total->Add(h2_i); outROOT->ReOpen("UPDATE"); outROOT->cd("accp_plots"); h2_total->Write("", TObject::kOverwrite); outROOT->ReOpen("READ");	       
	  
	}
	
      }//end loop over accp_HList

    //--------------------------------------------------------------------------------------------


    
    outROOT->Close();
    
  } //end else statement (assumes combined ROOTfile already exists)
  
}


//--------------------------MAIN ANALYSIS FUNCTIONS-----------------------------
void baseAnalyzer::run_data_analysis()
{
  /*
    Brief:  This method call all the necessary methods to carry out the 
    full data analysis. The main controls input parameters and cuts 
    are located in: set_basic_cuts.inp. The file to set the histogram
    bining is: set_basic_histos.inp. 

    This is supposed to be a generic baseAnalyzer class which analyzes data in
    a generic way. The analyzer assumes that all the calibrations from the data
    have been done. See each of the methods for details of what the analyzer does.

    Additional methods may be added in accordance with the necessity of the experiment.
    For example, methods to apply radiative corrections and get cross section still need 
    to be added, as well as a methods to determine the helicity data for beam asymmetry analysis. 

  */

  //------------------
  ReadInputFile();
  ReadReport();
  SetHistBins();
  CreateHist();
  
  ReadScalerTree();   
  ScalerEventLoop();       

  ReadTree();
  EventLoop();

  CalcEff();
  ApplyWeight();

  WriteHist();
  WriteReport();
  CombineHistos();

  //------------------
  
}
