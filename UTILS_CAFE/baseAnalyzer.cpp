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
baseAnalyzer::baseAnalyzer( int irun=-1, int ievt=-1, string mode="", string earm="", Bool_t ana_data=0, string ana_cuts="", string ana_type="", Bool_t hel_flag=0, string bcm_name="", double thrs=-1, string trig="", Bool_t combine_flag=0 )
  : run(irun), evtNum(ievt), daq_mode(mode), e_arm_name(earm), analyze_data(ana_data), analysis_cut(ana_cuts), analysis_type(ana_type), helicity_flag(hel_flag), bcm_type(bcm_name), bcm_thrs(thrs), trig_type(trig), combine_runs_flag(combine_flag)   //initialize member list 
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
  rand_HList = NULL;
  randSub_HList = NULL;
  
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
  H_ep_ctime_total  = NULL;
  H_ep_ctime  = NULL;
  
  //-HMS-
  H_hCerNpeSum      = NULL;  
  H_hCalEtotNorm    = NULL;
  H_hCalEtotTrkNorm = NULL;
  H_hHodBetaNtrk    = NULL;   
  H_hHodBetaTrk    = NULL;
		  
  //-SHMS-
  H_pNGCerNpeSum    = NULL;
  H_pHGCerNpeSum    = NULL;
  H_pCalEtotNorm    = NULL;
  H_pCalEtotTrkNorm = NULL;
  H_pHodBetaNtrk    = NULL;   
  H_pHodBetaTrk     = NULL;
  
  //HMS 2D PID
  H_hcal_vs_hcer     = NULL;  
  
  //sHMS 2D PID
  H_pcal_vs_phgcer   = NULL;  
  H_pcal_vs_pngcer   = NULL;  
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
  H_Em_src     = NULL;	
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

  // 2D Kinematics Histos
  H_Em_nuc_vs_Pm = NULL;
  H_Em_src_vs_Pm = NULL;
  
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


  //----------------------------------
  // Selected Histograms for Random
  // Coincidence Background Subtraction
  //----------------------------------

  H_ep_ctime_rand      =  NULL; 
  H_ep_ctime_rand_sub  =  NULL;     
                       		
  H_W_rand  =  NULL;	       	
  H_W_rand_sub  =  NULL;	
                       		       
  H_Q2_rand  =  NULL;	       	
  H_Q2_rand_sub  =  NULL;       
                       		
  H_xbj_rand  =  NULL;	       	
  H_xbj_rand_sub  =  NULL;      
                       		
  H_nu_rand  =  NULL;	       	
  H_nu_rand_sub  =  NULL;       
                       		
  H_q_rand  =  NULL;	       	
  H_q_rand_sub  =  NULL;	
                       		       
  H_Em_rand  =  NULL;	       	
  H_Em_rand_sub  =  NULL;       
                       		
  H_Em_nuc_rand  =  NULL;       
  H_Em_nuc_rand_sub  =  NULL;   
                       		
  H_Pm_rand  =  NULL;	       	
  H_Pm_rand_sub  =  NULL;       
                       		
  H_MM_rand  =  NULL;	       	
  H_MM_rand_sub  =  NULL;       
                       		
  H_thxq_rand  =  NULL;	       	
  H_thxq_rand_sub  =  NULL;     
                       		
  H_thrq_rand  =  NULL;	       	
  H_thrq_rand_sub = NULL;     	
  
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
  delete pid_HList;     pid_HList = NULL;
  delete kin_HList;     kin_HList = NULL;
  delete accp_HList;    accp_HList = NULL;
  delete rand_HList;    rand_HList = NULL;
  delete randSub_HList; randSub_HList = NULL;
  
  //-----------------------------
  // Detector Histogram Pointers
  //-----------------------------

  //-Coin. Time-
  delete H_ep_ctime_total; H_ep_ctime_total   = NULL;
  delete H_ep_ctime; H_ep_ctime   = NULL;

  //-HMS-
  delete H_hCerNpeSum;      H_hCerNpeSum      = NULL;  
  delete H_hCalEtotNorm;    H_hCalEtotNorm    = NULL;
  delete H_hCalEtotTrkNorm; H_hCalEtotTrkNorm = NULL;
  delete H_hHodBetaNtrk;    H_hHodBetaNtrk    = NULL;   
  delete H_hHodBetaTrk;     H_hHodBetaTrk     = NULL;
		  
  //-SHMS-
  delete H_pNGCerNpeSum;    H_pNGCerNpeSum    = NULL;
  delete H_pHGCerNpeSum;    H_pHGCerNpeSum    = NULL;
  delete H_pCalEtotNorm;    H_pCalEtotNorm    = NULL;
  delete H_pCalEtotTrkNorm; H_pCalEtotTrkNorm = NULL;
  delete H_pHodBetaNtrk;    H_pHodBetaNtrk    = NULL;   
  delete H_pHodBetaTrk;     H_pHodBetaTrk     = NULL;   
  
  //HMS 2D PID                      
  delete H_hcal_vs_hcer;         H_hcal_vs_hcer     = NULL;
  		     	     			    
  //sHMS 2D PID	     	     		    
  delete H_pcal_vs_phgcer;   	 H_pcal_vs_phgcer   = NULL;
  delete H_pcal_vs_pngcer;   	 H_pcal_vs_pngcer   = NULL;
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
  delete H_Em_src;	   H_Em_src     = NULL;	
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

  // 2D Kinematics Histos
  delete H_Em_nuc_vs_Pm;    H_Em_nuc_vs_Pm = NULL;
  delete H_Em_src_vs_Pm;    H_Em_src_vs_Pm = NULL;
  
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

  //----------------------------------
  // Selected Histograms for Random
  // Coincidence Background Subtraction
  //----------------------------------
  
  delete H_ep_ctime_rand;         H_ep_ctime_rand      =  NULL; 
  delete H_ep_ctime_rand_sub; 	  H_ep_ctime_rand_sub  =  NULL; 
                       		                       		
  delete H_W_rand;	       	  H_W_rand  =  NULL;	       	
  delete H_W_rand_sub;	       	  H_W_rand_sub  =  NULL;	
                       		                       		
  delete H_Q2_rand;	       	  H_Q2_rand  =  NULL;	       	
  delete H_Q2_rand_sub;       	  H_Q2_rand_sub  =  NULL;       
                       		                       		
  delete H_xbj_rand;	       	  H_xbj_rand  =  NULL;	       	
  delete H_xbj_rand_sub;      	  H_xbj_rand_sub  =  NULL;      
                       		                       		
  delete H_nu_rand;	       	  H_nu_rand  =  NULL;	       	
  delete H_nu_rand_sub;       	  H_nu_rand_sub  =  NULL;       
                       		                       		
  delete H_q_rand;	       	  H_q_rand  =  NULL;	       	
  delete H_q_rand_sub;	       	  H_q_rand_sub  =  NULL;	
                       		                       		
  delete H_Em_rand;	       	  H_Em_rand  =  NULL;	       	
  delete H_Em_rand_sub;       	  H_Em_rand_sub  =  NULL;       
                       		                       		
  delete H_Em_nuc_rand;       	  H_Em_nuc_rand  =  NULL;       
  delete H_Em_nuc_rand_sub;   	  H_Em_nuc_rand_sub  =  NULL;   
                       		                       		
  delete H_Pm_rand;	       	  H_Pm_rand  =  NULL;	       	
  delete H_Pm_rand_sub;       	  H_Pm_rand_sub  =  NULL;       
                       		                       		
  delete H_MM_rand;	       	  H_MM_rand  =  NULL;	       	
  delete H_MM_rand_sub;       	  H_MM_rand_sub  =  NULL;       
                       		                       		
  delete H_thxq_rand;	       	  H_thxq_rand  =  NULL;	       	
  delete H_thxq_rand_sub;     	  H_thxq_rand_sub  =  NULL;     
                       		                       		
  delete H_thrq_rand;	       	  H_thrq_rand  =  NULL;	       	
  delete H_thrq_rand_sub;     	  H_thrq_rand_sub = NULL;     	


  
  
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

  input_FileNamePattern = "UTILS_CAFE/inp/set_basic_filenames.inp";
  input_CutFileName     = "UTILS_CAFE/inp/set_basic_cuts.inp";
  input_HBinFileName    = "UTILS_CAFE/inp/set_basic_histos.inp";
  input_SIMCinfo_FileName = "UTILS_CAFE/inp/set_basic_simc_param.inp";
  //==========================================
  //     READ FILE NAME PATTERN
  //==========================================

  TString temp; //temporary string placeholder


  
  //-------------------------------
  //----INPUTS (USER READS IN)-----
  //-------------------------------
  
  //Define Input (.root) File Name Patterns (read principal ROOTfile from experiment)
  temp = trim(split(FindString("input_ROOTfilePattern", input_FileNamePattern.Data())[0], '=')[1]);
  data_InputFileName = Form(temp.Data(),  analysis_type.Data(), analysis_type.Data(), run, evtNum);

  //Check if ROOTfile exists
  in_file.open(data_InputFileName.Data());
  cout << "in_file.fail() --> " << in_file.fail() << endl;
  if(in_file.fail()){
    cout << Form("ROOTFile: %s does NOT exist ! ! !", data_InputFileName.Data()) << endl;
    cout << "Exiting NOW !" << endl;
    gSystem->Exit(0);
  }  
  in_file.close();
  
  //Define Input (.report) File Name Pattern (read principal REPORTfile from experiment)
  temp = trim(split(FindString("input_REPORTPattern", input_FileNamePattern.Data())[0], '=')[1]);
  data_InputReport = Form(temp.Data(), analysis_type.Data(), analysis_type.Data(), run, evtNum);
  
  //Check if REPORTFile exists
  in_file.open(data_InputReport.Data());
  cout << "in_file.fail() --> " << in_file.fail() << endl;
  if(in_file.fail()){
    cout << Form("REPORTFile: %s does NOT exist ! ! !", data_InputReport.Data()) << endl;
    cout << "Exiting NOW !" << endl;
    gSystem->Exit(0);
  }
  in_file.close();
  
  //----------------------------------
  //----OUTPUTS (USER WRITES OUT)-----
  //----------------------------------
  //Define Output (.root) File Name Pattern (analyzed histos are written to this file)
  temp = trim(split(FindString("output_ROOTfilePattern", input_FileNamePattern.Data())[0], '=')[1]);
  data_OutputFileName = Form(temp.Data(), analysis_type.Data(), run, evtNum);

  //Define Output (.root) File Name Pattern (analyzed combined histos are written to this file)
  temp = trim(split(FindString("output_ROOTfilePattern_final", input_FileNamePattern.Data())[0], '=')[1]);
  data_OutputFileName_combined = Form(temp.Data(), analysis_type.Data());

  //Define Output (.txt) File Name Pattern (analysis report is written to this file) -- append run numbers
  temp = trim(split(FindString("output_SummaryPattern", input_FileNamePattern.Data())[0], '=')[1]);
  output_SummaryFileName = Form(temp.Data(), analysis_type.Data());

  //Define Output (.txt) File Name Pattern (analysis report is written to this file) -- short report on a per-run basis
  temp = trim(split(FindString("output_REPORTPattern", input_FileNamePattern.Data())[0], '=')[1]);
  output_ReportFileName = Form(temp.Data(), analysis_type.Data(), run, evtNum);


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
  
  ePctime_cut_flag = stoi(split(FindString("ePctime_cut_flag", input_CutFileName.Data())[0], '=')[1]);

  // main coincidence time peak min/max window cut
  ePctime_cut_min = stod(split(FindString("ePctime_cut_min", input_CutFileName.Data())[0], '=')[1]);
  ePctime_cut_max = stod(split(FindString("ePctime_cut_max", input_CutFileName.Data())[0], '=')[1]);

  // accidentals to the right of main coin. peak
  ePctime_cut_min_R = stod(split(FindString("ePctime_cut_min_R", input_CutFileName.Data())[0], '=')[1]);
  ePctime_cut_max_R = stod(split(FindString("ePctime_cut_max_R", input_CutFileName.Data())[0], '=')[1]);

  // accidentals to the left of main coin. peak
  ePctime_cut_min_L = stod(split(FindString("ePctime_cut_min_L", input_CutFileName.Data())[0], '=')[1]);
  ePctime_cut_max_L = stod(split(FindString("ePctime_cut_max_L", input_CutFileName.Data())[0], '=')[1]);

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
  

  //(HMS PID) Calorimeter Total Energy Normalized By Track Momentum
  hetot_trkNorm_pidCut_flag = stoi(split(FindString("hetot_trkNorm_pidCut_flag", input_CutFileName.Data())[0], '=')[1]);
  cpid_hetot_trkNorm_min = stod(split(FindString("cpid_hetot_trkNorm_min", input_CutFileName.Data())[0], '=')[1]);
  cpid_hetot_trkNorm_max = stod(split(FindString("cpid_hetot_trkNorm_max", input_CutFileName.Data())[0], '=')[1]);
  
  //(HMS PID) Gas Cherenkov
  hcer_pidCut_flag = stoi(split(FindString("hcer_pidCut_flag", input_CutFileName.Data())[0], '=')[1]);
  cpid_hcer_npeSum_min = stod(split(FindString("cpid_hcer_npeSum_min", input_CutFileName.Data())[0], '=')[1]);
  cpid_hcer_npeSum_max = stod(split(FindString("cpid_hcer_npeSum_max", input_CutFileName.Data())[0], '=')[1]);


  //-----Kinematics Cuts------
  // H(e,e'p)
  
  //4-Momentum Transfers [GeV^2]
  Q2_heep_cut_flag = stoi(split(FindString("Q2_heep_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_heep_Q2_min = stod(split(FindString("c_heep_Q2_min", input_CutFileName.Data())[0], '=')[1]);
  c_heep_Q2_max = stod(split(FindString("c_heep_Q2_max", input_CutFileName.Data())[0], '=')[1]);

  //bjorken-x
  xbj_heep_cut_flag = stoi(split(FindString("xbj_heep_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_heep_xbj_min = stod(split(FindString("c_heep_xbj_min", input_CutFileName.Data())[0], '=')[1]);
  c_heep_xbj_max = stod(split(FindString("c_heep_xbj_max", input_CutFileName.Data())[0], '=')[1]);

  
  //Missing Energy [GeV]
  Em_heep_cut_flag = stoi(split(FindString("Em_heep_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_heep_Em_min = stod(split(FindString("c_heep_Em_min", input_CutFileName.Data())[0], '=')[1]);
  c_heep_Em_max = stod(split(FindString("c_heep_Em_max", input_CutFileName.Data())[0], '=')[1]);
  
  //Invariant Mass, W [GeV]
  W_heep_cut_flag = stoi(split(FindString("W_heep_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_heep_W_min = stod(split(FindString("c_heep_W_min", input_CutFileName.Data())[0], '=')[1]);
  c_heep_W_max = stod(split(FindString("c_heep_W_max", input_CutFileName.Data())[0], '=')[1]);
  
  //Missing Mass Cut (Check which MM Cut is actually being applied: By default, it should be proton MM) 
  //Protons
  MM_heep_cut_flag = stoi(split(FindString("MM_heep_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_heep_MM_min = stod(split(FindString("c_heep_MM_min", input_CutFileName.Data())[0], '=')[1]);
  c_heep_MM_max = stod(split(FindString("c_heep_MM_max", input_CutFileName.Data())[0], '=')[1]);

  // CaFe A(e,e'p) Mean-Field (MF) Kinematic Cuts
  // 4-Momentum Transfers [GeV^2]
  Q2_MF_cut_flag = stoi(split(FindString("Q2_MF_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_MF_Q2_min = stod(split(FindString("c_MF_Q2_min", input_CutFileName.Data())[0], '=')[1]);
  c_MF_Q2_max = stod(split(FindString("c_MF_Q2_max", input_CutFileName.Data())[0], '=')[1]);
  
  // Missing Momentum [GeV]
  Pm_MF_cut_flag = stoi(split(FindString("Pm_MF_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_MF_Pm_min = stod(split(FindString("c_MF_Pm_min", input_CutFileName.Data())[0], '=')[1]);
  c_MF_Pm_max = stod(split(FindString("c_MF_Pm_max", input_CutFileName.Data())[0], '=')[1]);

  // Missing Energy [GeV] --- ONLY for deuteron target
  Em_d2MF_cut_flag = stoi(split(FindString("Em_d2MF_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_d2MF_Em_min = stod(split(FindString("c_d2MF_Em_min", input_CutFileName.Data())[0], '=')[1]);
  c_d2MF_Em_max = stod(split(FindString("c_d2MF_Em_max", input_CutFileName.Data())[0], '=')[1]);

  // Missing Energy [GeV] --- ONLY for MF A>2 nuclei
  Em_MF_cut_flag = stoi(split(FindString("Em_MF_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_MF_Em_min = stod(split(FindString("c_MF_Em_min", input_CutFileName.Data())[0], '=')[1]);
  c_MF_Em_max = stod(split(FindString("c_MF_Em_max", input_CutFileName.Data())[0], '=')[1]);
  
  // CaFe A(e,e'p) Short-Range Correlations (SRC) Kinematic Cuts 
  // 4-Momentum Transfers [GeV^2]
  Q2_SRC_cut_flag = stoi(split(FindString("Q2_SRC_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_SRC_Q2_min = stod(split(FindString("c_SRC_Q2_min", input_CutFileName.Data())[0], '=')[1]);
  c_SRC_Q2_max = stod(split(FindString("c_SRC_Q2_max", input_CutFileName.Data())[0], '=')[1]);

  // Missing Momentum [GeV]
  Pm_SRC_cut_flag = stoi(split(FindString("Pm_SRC_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_SRC_Pm_min = stod(split(FindString("c_SRC_Pm_min", input_CutFileName.Data())[0], '=')[1]);
  c_SRC_Pm_max = stod(split(FindString("c_SRC_Pm_max", input_CutFileName.Data())[0], '=')[1]);

  // x-Bjorken
  Xbj_SRC_cut_flag = stoi(split(FindString("Xbj_SRC_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_SRC_Xbj_min = stod(split(FindString("c_SRC_Xbj_min", input_CutFileName.Data())[0], '=')[1]);
  c_SRC_Xbj_max = stod(split(FindString("c_SRC_Xbj_max", input_CutFileName.Data())[0], '=')[1]);

  // in-plane recoil (undetected) angle, theta_rq [deg]
  thrq_SRC_cut_flag = stoi(split(FindString("thrq_SRC_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_SRC_thrq_min = stod(split(FindString("c_SRC_thrq_min", input_CutFileName.Data())[0], '=')[1]);
  c_SRC_thrq_max = stod(split(FindString("c_SRC_thrq_max", input_CutFileName.Data())[0], '=')[1]);

  // Missing Energy [GeV] --- ONLY for deuteron target
  Em_d2SRC_cut_flag = stoi(split(FindString("Em_d2SRC_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_d2SRC_Em_min = stod(split(FindString("c_d2SRC_Em_min", input_CutFileName.Data())[0], '=')[1]);
  c_d2SRC_Em_max = stod(split(FindString("c_d2SRC_Em_max", input_CutFileName.Data())[0], '=')[1]);

  // Missing Energy [GeV] --- ONLY for MF A>2 nuclei
  Em_SRC_cut_flag = stoi(split(FindString("Em_SRC_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  
  //------Acceptance Cuts-------
  
  //Hadron Arm
  hdelta_cut_flag = stoi(split(FindString("hdelta_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_hdelta_min = stod(split(FindString("c_hdelta_min", input_CutFileName.Data())[0], '=')[1]);
  c_hdelta_max = stod(split(FindString("c_hdelta_max", input_CutFileName.Data())[0], '=')[1]);

  hxptar_cut_flag = stoi(split(FindString("hxptar_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_hxptar_min = stod(split(FindString("c_hxptar_min", input_CutFileName.Data())[0], '=')[1]);
  c_hxptar_max = stod(split(FindString("c_hxptar_max", input_CutFileName.Data())[0], '=')[1]);

  hyptar_cut_flag = stoi(split(FindString("hyptar_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_hyptar_min = stod(split(FindString("c_hyptar_min", input_CutFileName.Data())[0], '=')[1]);
  c_hyptar_max = stod(split(FindString("c_hyptar_max", input_CutFileName.Data())[0], '=')[1]);

  hmsCollCut_flag = stoi(split(FindString("hmsCollCut_flag", input_CutFileName.Data())[0], '=')[1]);
  
  //Electron Arm
  edelta_cut_flag = stoi(split(FindString("edelta_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_edelta_min = stod(split(FindString("c_edelta_min", input_CutFileName.Data())[0], '=')[1]);
  c_edelta_max = stod(split(FindString("c_edelta_max", input_CutFileName.Data())[0], '=')[1]);

  exptar_cut_flag = stoi(split(FindString("exptar_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_exptar_min = stod(split(FindString("c_exptar_min", input_CutFileName.Data())[0], '=')[1]);
  c_exptar_max = stod(split(FindString("c_exptar_max", input_CutFileName.Data())[0], '=')[1]);

  eyptar_cut_flag = stoi(split(FindString("eyptar_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_eyptar_min = stod(split(FindString("c_eyptar_min", input_CutFileName.Data())[0], '=')[1]);
  c_eyptar_max = stod(split(FindString("c_eyptar_max", input_CutFileName.Data())[0], '=')[1]);

  shmsCollCut_flag = stoi(split(FindString("shmsCollCut_flag", input_CutFileName.Data())[0], '=')[1]);

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
  double temp_var;

  // Read Target Mass
  temp = FindString("Target_Mass_amu",  data_InputReport.Data())[0];
  temp_var = stod(split(temp, ':')[1]);

  double max_diff = 1e-6;
  
  if(abs(temp_var-MH_amu)<=max_diff){
    tgt_type = "LH2";
    tgt_mass = MH_amu; tgt_density = rho_H; tgt_thickness = thick_H;

    
  }
  
  else if(abs(temp_var-MD_amu)<=max_diff){
    tgt_type = "LD2";
    tgt_mass = MD_amu; tgt_density = rho_D; tgt_thickness = thick_D;

  }
  
  else if(abs(temp_var-MBe9_amu)<=max_diff){
    tgt_type = "Be9";
    tgt_mass = MBe9_amu; tgt_density = rho_Be9; tgt_thickness = thick_Be9;    
  }
  
  else if(abs(temp_var-MB10_amu)<=max_diff){
    tgt_type = "B10";
    tgt_mass = MB10_amu; tgt_density = rho_B10; tgt_thickness = thick_B10;
  }
  
  else if(abs(temp_var-MB11_amu)<=max_diff){
    tgt_type = "B11";
    tgt_mass = MB11_amu; tgt_density = rho_B11; tgt_thickness = thick_B11;
  }
 
  else if(abs(temp_var-MC12_amu)<=max_diff){    
    if(analysis_cut=="optics") {tgt_type = "C12_optics";}
    else{tgt_type = "C12";}
    tgt_mass = MC12_amu; tgt_density = rho_C12; tgt_thickness = thick_C12;
  }
  
  else if(abs(temp_var-MAl27_amu)<=max_diff){ 
    tgt_type = "Al27";
    tgt_mass = MAl27_amu;  tgt_density = rho_Al27; tgt_thickness = thick_Al27;
  }

  else if(abs(temp_var-MCa40_amu)<=max_diff){
    tgt_type = "Ca40";
    tgt_mass = MCa40_amu; tgt_density = rho_Ca40; tgt_thickness = thick_Ca40;
  }

  else if(abs(temp_var-MCa48_amu)<=max_diff){
    tgt_type = "Ca48";
    tgt_mass = MCa48_amu; tgt_density = rho_Ca48; tgt_thickness = thick_Ca48;
  }

  else if(abs(temp_var-MFe54_amu)<=max_diff){
    tgt_type = "Fe54"; 
    tgt_mass = MFe54_amu; tgt_density = rho_Fe54; tgt_thickness = thick_Fe54;
  }

  else if(abs(temp_var-MTi48_amu)<=max_diff){
    tgt_type = "Ti48";
    tgt_mass = MTi48_amu; tgt_density = rho_Ti48; tgt_thickness = thick_Ti48;
  }
  
  else{
    cout << "Target mass (amu) mis-match of >1E-6 between this script and standard.kinematics . . . Check target mass is set correctly in standard.kinematics file !" << endl;
    gSystem->Exit(0);
    
  }
  
  //Read Pre-Scale Factors
  temp =  FindString("Ps1_factor", data_InputReport.Data())[0]; 
  Ps1_factor = stod(split(temp, ':')[1]);

  temp =  FindString("Ps2_factor", data_InputReport.Data())[0]; 
  Ps2_factor = stod(split(temp, ':')[1]);

  temp =  FindString("Ps3_factor", data_InputReport.Data())[0]; 
  Ps3_factor = stod(split(temp, ':')[1]);
  
  temp =  FindString("Ps4_factor", data_InputReport.Data())[0]; 
  Ps4_factor = stod(split(temp, ':')[1]);
 
  temp =  FindString("Ps5_factor", data_InputReport.Data())[0]; 
  Ps5_factor = stod(split(temp, ':')[1]);
  
  temp =  FindString("Ps6_factor", data_InputReport.Data())[0]; 
  Ps6_factor = stod(split(temp, ':')[1]);
  

  //Read spec. kinematics
  temp = FindString("Beam_Energy",  data_InputReport.Data())[0];
  beam_energy = stod(split(temp, ':')[1]);

  //hms spec. kinematics
  temp = FindString("HMS_Particle_Mass",  data_InputReport.Data())[0];  // GeV
  hms_part_mass = stod(split(temp, ':')[1]);

  temp = FindString("HMS_P_Central",  data_InputReport.Data())[0]; //GeV/v
  hms_p = stod(split(temp, ':')[1]);
  
  temp = FindString("HMS_Angle",  data_InputReport.Data())[0];  // deg
  hms_angle = stod(split(temp, ':')[1]);

  // shms spec. kinematics
  temp = FindString("SHMS_Particle_Mass",  data_InputReport.Data())[0];
  shms_part_mass = stod(split(temp, ':')[1]);

  temp = FindString("SHMS_P_Central",  data_InputReport.Data())[0];
  shms_p = stod(split(temp, ':')[1]);
  
  temp = FindString("SHMS_Angle",  data_InputReport.Data())[0];
  shms_angle = stod(split(temp, ':')[1]);

  // run start_time (format: yyyy-mm-dd HH:MM:SS)
  temp = FindString("start_of_run", data_InputReport.Data())[0];
  start_of_run =  split(temp, '=')[1];

  // run end_time (format: yyyy-mm-dd HH:MM:SS)
  temp = FindString("end_of_run", data_InputReport.Data())[0];
  end_of_run =  split(temp, '=')[1];
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

  Em_nuc_nbins     	 = stod(split(FindString("Em_nuc_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  Em_nuc_xmin      	 = stod(split(FindString("Em_nuc_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  Em_nuc_xmax      	 = stod(split(FindString("Em_nuc_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  
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

  rand_HList = new TList();
  randSub_HList = new TList();
  
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
  H_ep_ctime   = new TH1F("H_ep_ctime", "ep Coincidence Time; ep Coincidence Time [ns]; Counts ", coin_nbins, coin_xmin, coin_xmax);
  H_ep_ctime->Sumw2(); //Apply sum of weight squared to this histogram ABOVE.
  H_ep_ctime->SetDefaultSumw2(kTRUE);  //Generalize sum weights squared to all histograms  (ROOT 6 has this by default. ROOT 5 does NOT)
  H_ep_ctime_total   = new TH1F("H_ep_ctime_total", "ep Coincidence Time; ep Coincidence Time [ns]; Counts ", coin_nbins, coin_xmin, coin_xmax);


  //HMS DETECTORS HISTOS
  H_hCerNpeSum      = new TH1F("H_hCerNpeSum", "HMS Cherenkov NPE Sum; Cherenkov NPE Sum; Counts ", hcer_nbins, hcer_xmin, hcer_xmax);
  H_hCalEtotNorm    = new TH1F("H_hCalEtotNorm", "HMS Calorimeter Normalized Total Energy; E_{tot} / P_{cent}; Counts ", hcal_nbins, hcal_xmin, hcal_xmax);
  H_hCalEtotTrkNorm = new TH1F("H_hCalEtotTrkNorm", "HMS Calorimeter Total Normalized Track Energy; E_{tot} / P_{trk}; Counts ", hcal_nbins, hcal_xmin, hcal_xmax);
  H_hHodBetaNtrk    = new TH1F("H_hHodBetaNtrk", "HMS Hodo #beta (no track); #beta (no track); Counts ", hbeta_nbins, hbeta_xmin, hbeta_xmax);
  H_hHodBetaTrk     = new TH1F("H_hHodBetaTrk", "HMS Hodo #beta (golden track); #beta (golden track); Counts ", hbeta_nbins, hbeta_xmin, hbeta_xmax);
  
  //SHMS DETECTORS HISTOS
  H_pNGCerNpeSum    = new TH1F("H_pNGCerNpeSum", "SHMS Noble Gas Cherenkov NPE Sum; Cherenkov NPE Sum; Counts  ", pngcer_nbins, pngcer_xmin, pngcer_xmax);
  H_pHGCerNpeSum    = new TH1F("H_pHGCerNpeSum", "SHMS Heavy Gas Cherenkov NPE Sum; Cherenkov NPE Sum; Counts  ", phgcer_nbins, phgcer_xmin, phgcer_xmax);
  H_pCalEtotNorm    = new TH1F("H_pCalEtotNorm", "SHMS Calorimeter Normalized Total Energy; E_{tot} / P_{cent}; Counts ", pcal_nbins, pcal_xmin, pcal_xmax);
  H_pCalEtotTrkNorm = new TH1F("H_pCalEtotTrkNorm", "SHMS Calorimeter Total Normalized Track Energy; E_{tot} / P_{trk}; Counts ", pcal_nbins, pcal_xmin, pcal_xmax);
  H_pHodBetaNtrk    = new TH1F("H_pBetaNtrk", "SHMS Hodo #beta (no track); #beta (no track); Counts ", pbeta_nbins, pbeta_xmin, pbeta_xmax);
  H_pHodBetaTrk     = new TH1F("H_pBetaTrk", "SHMS Hodo #beta (golden track); #beta (golden track); Counts ", pbeta_nbins, pbeta_xmin, pbeta_xmax);

  //HMS 2D PID               
  H_hcal_vs_hcer     = new TH2F("H_hcal_vs_hcer", "HMS: Calorimeter vs. Cherenkov; Calorimeter E_{tot}/P_{trk}; Cherenkov NPE Sum", hcal_nbins, hcal_xmin, hcal_xmax, hcer_nbins, hcer_xmin, hcer_xmax);     
  		     	     
  //SHMS 2D PID	     	     
  H_pcal_vs_phgcer    = new TH2F("H_pcal_vs_phgcer", "SHMS: Heavy Gas Cherenkov (HGC) vs. Calorimeter; Calorimeter E_{tot}/P_{trk}; HGC NPE Sum", pcal_nbins, pcal_xmin, pcal_xmax, phgcer_nbins, phgcer_xmin, phgcer_xmax);        
  H_pcal_vs_pngcer    = new TH2F("H_pcal_vs_pngcer", "SHMS: Noble Gas Cherenkov (NGC) vs. Calorimeter; Calorimeter E_{tot}/P_{trk}; NGC NPE Sum", pcal_nbins, pcal_xmin, pcal_xmax, pngcer_nbins, pngcer_xmin, pngcer_xmax);   
  H_pngcer_vs_phgcer  = new TH2F("H_pngcer_vs_phgcer","SHMS: Heavy Gas Cherenkov (HGC) vs. Noble Gas Cherenkov (NGC); NGC NPE Sum; HGC NPE Sum", pngcer_nbins, pngcer_xmin, pngcer_xmax, phgcer_nbins, phgcer_xmin, phgcer_xmax); 
  
  
  //Add PID Histos to TList
  pid_HList->Add(H_ep_ctime_total);
  pid_HList->Add(H_ep_ctime);
  pid_HList->Add(H_hCerNpeSum);
  pid_HList->Add(H_hCalEtotNorm);
  pid_HList->Add(H_hCalEtotTrkNorm);
  pid_HList->Add(H_hHodBetaNtrk);
  pid_HList->Add(H_hHodBetaTrk);
  pid_HList->Add(H_pNGCerNpeSum);
  pid_HList->Add(H_pHGCerNpeSum);
  pid_HList->Add(H_pCalEtotNorm);
  pid_HList->Add(H_pCalEtotTrkNorm);
  pid_HList->Add(H_pHodBetaNtrk);
  pid_HList->Add(H_pHodBetaTrk);
  //Add 2D PID
  pid_HList->Add(H_hcal_vs_hcer);  
  pid_HList->Add(H_pcal_vs_phgcer);
  pid_HList->Add(H_pcal_vs_pngcer);
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
  H_thq     = new TH1F("H_thq", "#theta_{q}", thq_nbins, thq_xmin, thq_xmax); 
  H_phq     = new TH1F("H_phq", "#phi_{q}", phq_nbins, phq_xmin, phq_xmax); 
  H_epsilon = new TH1F("H_epsilon", "Virtual Photon (#gamma) Polarization Factor" , epsilon_nbins, epsilon_xmin, epsilon_xmax); 

  //Secondary (Hadron) Kinematics (recoil and missing are used interchageably) ()
  H_Em      = new TH1F("H_Em","Missing Energy", Em_nbins, Em_xmin, Em_xmax);   
  H_Em_nuc  = new TH1F("H_Em_nuc","Nuclear Missing Energy", Em_nuc_nbins, Em_nuc_xmin, Em_nuc_xmax);
  H_Em_src  = new TH1F("H_Em_src","SRC Nuclear Missing Energy", Em_nuc_nbins, Em_nuc_xmin, Em_nuc_xmax); 
  H_Pm      = new TH1F("H_Pm","Missing Momentum, P_{miss}", Pm_nbins, Pm_xmin, Pm_xmax); 
  H_Pmx_lab = new TH1F("H_Pmx_Lab","P_{miss, x} (Lab)", Pmx_lab_nbins, Pmx_lab_xmin, Pmx_lab_xmax);         
  H_Pmy_lab = new TH1F("H_Pmy_Lab","P_{miss, y} (Lab)", Pmy_lab_nbins, Pmy_lab_xmin, Pmy_lab_xmax);    
  H_Pmz_lab = new TH1F("H_Pmz_Lab","P_{miss, z} (Lab)", Pmz_lab_nbins, Pmz_lab_xmin, Pmz_lab_xmax);  
  H_Pmx_q   = new TH1F("H_Pmx_q","P_{miss, xq} (w.r.t #vec{q}) ", Pmx_q_nbins, Pmx_q_xmin, Pmx_q_xmax);   
  H_Pmy_q   = new TH1F("H_Pmy_q","P_{miss, yq} (w.r.t #vec{q}) ", Pmy_q_nbins, Pmy_q_xmin, Pmy_q_xmax); 
  H_Pmz_q   = new TH1F("H_Pmz_q","P_{miss, zq} (along #vec{q}) ", Pmz_q_nbins, Pmz_q_xmin, Pmz_q_xmax); 
  H_Tx      = new TH1F("H_Tx", "Kinetic Energy, T_{p} (detected)", Tx_nbins, Tx_xmin, Tx_xmax);     
  H_Tr      = new TH1F("H_Tr", "Kinetic Energy, T_{r} (recoil)",   Tr_nbins, Tr_xmin, Tr_xmax);  
  H_MM      = new TH1F("H_MM","Missing Mass, M_{miss}", MM_nbins, MM_xmin, MM_xmax);        
  H_MM2     = new TH1F("H_MM2","Missing Mass Squared, M^{2}_{miss}", MM2_nbins, MM2_xmin, MM2_xmax); 
  H_thx     = new TH1F("H_thx", "Hadron Scattering Angle (detected), #theta_{p}", thx_nbins, thx_xmin, thx_xmax);
  H_Pf      = new TH1F("H_Pf", "Final Hadron Momentum (detected), p_{f}", Pf_nbins, Pf_xmin, Pf_xmax);
  H_thxq    = new TH1F("H_thxq", "In-Plane (detected) Angle, #theta_{pq}", thxq_nbins, thxq_xmin, thxq_xmax);
  H_thrq    = new TH1F("H_thrq", "In-Plane (recoil) Angle, #theta_{rq}", thrq_nbins, thrq_xmin, thrq_xmax);
  H_phxq    = new TH1F("H_phxq", "Out-of-Plane (detected) Angle, #phi_{pq}", phxq_nbins, phxq_xmin, phxq_xmax);
  H_phrq    = new TH1F("H_phrq", "Out-of-Plane (recoil) Angle, #phi_{rq}", phrq_nbins, phrq_xmin, phrq_xmax);
  H_Tx_cm   = new TH1F("H_Tx_cm", "Kinetic Energy, T_{x, cm} (detected)", Tx_cm_nbins, Tx_cm_xmin, Tx_cm_xmax);     
  H_Tr_cm   = new TH1F("H_Tr_cm", "Kinetic Energy, T_{r, cm} (recoil)",   Tr_cm_nbins, Tr_cm_xmin, Tr_cm_xmax);  
  H_thxq_cm = new TH1F("H_thxq_cm", "In-Plane (detected) Angle, #theta_{pq, cm}", thxq_cm_nbins, thxq_cm_xmin, thxq_cm_xmax);
  H_thrq_cm = new TH1F("H_thrq_cm", "In-Plane (recoil) Angle, #theta_{rq, cm}", thrq_cm_nbins, thrq_cm_xmin, thrq_cm_xmax);
  H_phxq_cm = new TH1F("H_phxq_cm", "Out-of-Plane (detected) Angle, #phi_{pq, cm}", phxq_cm_nbins, phxq_cm_xmin, phxq_cm_xmax);
  H_phrq_cm = new TH1F("H_phrq_cm", "Out-of-Plane (recoil) Angle, #phi_{rq, cm}", phrq_cm_nbins, phrq_cm_xmin, phrq_cm_xmax);
  H_Ttot_cm = new TH1F("H_Ttot_cm", "Total CM Kinetic Energy, T_{tot,cm}", Ttot_cm_nbins, Ttot_cm_xmin, Ttot_cm_xmax);
  H_MandelS = new TH1F("H_MandelS", "s-Mandelstam", MandelS_nbins, MandelS_xmin, MandelS_xmax);     
  H_MandelT = new TH1F("H_MandelT", "t-Mandelstam", MandelT_nbins, MandelT_xmin, MandelT_xmax);
  H_MandelU = new TH1F("H_MandelU", "u-Mandelstam", MandelU_nbins, MandelU_xmin, MandelU_xmax);     
  
  // (Cosine, Sine) Histos of detected AND recoil angles (range is fixed at: -1, 1)
  //LAB FRAME
  H_cth_xq = new TH1F("H_cth_xq", "cos(#theta_{pq})", thxq_nbins, -1, 1);
  H_cth_rq = new TH1F("H_cth_rq", "cos(#theta_{rq})", thrq_nbins, -1, 1);
  H_sth_xq = new TH1F("H_sth_xq", "sin(#theta_{pq})", thxq_nbins, -1, 1);
  H_sth_rq = new TH1F("H_sth_rq", "sin(#theta_{rq})", thrq_nbins, -1, 1);
  H_cphi_xq = new TH1F("H_cphi_xq", "cos(#phi_{pq})", phxq_nbins, -1, 1);
  H_cphi_rq = new TH1F("H_cphi_rq", "cos(#phi_{rq})", phrq_nbins, -1, 1);
  H_sphi_xq = new TH1F("H_sphi_xq", "sin(#phi_{pq})", phxq_nbins, -1, 1);
  H_sphi_rq = new TH1F("H_sphi_rq", "sin(#phi_{rq})", phrq_nbins, -1, 1);
  //CM FRAME
  H_cth_xq_cm = new TH1F("H_cth_xq_cm", "cos(#theta_{pq,cm})", thxq_cm_nbins, -1, 1);
  H_cth_rq_cm = new TH1F("H_cth_rq_cm", "cos(#theta_{rq,cm})", thrq_cm_nbins, -1, 1);
  H_sth_xq_cm = new TH1F("H_sth_xq_cm", "sin(#theta_{pq,cm})", thxq_cm_nbins, -1, 1);
  H_sth_rq_cm = new TH1F("H_sth_rq_cm", "sin(#theta_{rq,cm})", thrq_cm_nbins, -1, 1);
  H_cphi_xq_cm = new TH1F("H_cphi_xq_cm", "cos(#phi_{pq,cm})", phxq_cm_nbins, -1, 1);
  H_cphi_rq_cm = new TH1F("H_cphi_rq_cm", "cos(#phi_{rq,cm})", phrq_cm_nbins, -1, 1);
  H_sphi_xq_cm = new TH1F("H_sphi_xq_cm", "sin(#phi_{pq,cm})", phxq_cm_nbins, -1, 1);
  H_sphi_rq_cm = new TH1F("H_sphi_rq_cm", "sin(#phi_{rq,cm})", phrq_cm_nbins, -1, 1);

  // 2d kin histos
  H_Em_nuc_vs_Pm = new TH2F("H_Em_nuc_vs_Pm", "Em_nuc vs. Pm", Pm_nbins, Pm_xmin, Pm_xmax, Em_nuc_nbins, Em_nuc_xmin, Em_nuc_xmax);
  H_Em_src_vs_Pm = new TH2F("H_Em_src_vs_Pm", "Em_src vs. Pm", Pm_nbins, Pm_xmin, Pm_xmax, Em_nuc_nbins, Em_nuc_xmin, Em_nuc_xmax);
  
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
  kin_HList->Add( H_Em_src   );
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

  // 2d kin histos
  kin_HList->Add( H_Em_nuc_vs_Pm );
  kin_HList->Add( H_Em_src_vs_Pm );
  
  //----------------------------------------------------------------------
  //---------HISTOGRAM CATEGORY: Spectrometer Acceptance  (ACCP)----------
  //----------------------------------------------------------------------


  //Electron Arm Focal Plane Quantities
  H_exfp = new TH1F("H_exfp", Form("%s X_{fp}; X_{fp} [cm]; Counts ", e_arm_name.Data()), exfp_nbins, exfp_xmin, exfp_xmax);
  H_eyfp = new TH1F("H_eyfp", Form("%s Y_{fp}; Y_{fp} [cm]; Counts ", e_arm_name.Data()), eyfp_nbins, eyfp_xmin, eyfp_xmax);
  H_expfp = new TH1F("H_expfp", Form("%s X'_{fp}; X'_{fp} [rad]; Counts ", e_arm_name.Data()), expfp_nbins, expfp_xmin, expfp_xmax);
  H_eypfp = new TH1F("H_eypfp", Form("%s Y'_{fp}; Y'_{fp} [rad]; Counts ", e_arm_name.Data()), eypfp_nbins, eypfp_xmin, eypfp_xmax);
  
  //Electron Arm Reconstructed Quantities 
  H_eytar = new TH1F("H_eytar", Form("%s Y_{tar}; Y_{tar} [cm]; Counts ", e_arm_name.Data()), eytar_nbins, eytar_xmin, eytar_xmax);
  H_exptar = new TH1F("H_exptar", Form("%s X'_{tar}; X'_{tar} [rad]; Counts ", e_arm_name.Data()), exptar_nbins, exptar_xmin, exptar_xmax);
  H_eyptar = new TH1F("H_eyptar", Form("%s Y'_{tar}; Y'_{tar} [rad]; Counts ", e_arm_name.Data()), eyptar_nbins, eyptar_xmin, eyptar_xmax);
  H_edelta = new TH1F("H_edelta", Form("%s Momentum Acceptance, #delta; #delta [%%]; Counts ", e_arm_name.Data()), edelta_nbins, edelta_xmin, edelta_xmax);
  
  //Hadron arm Focal Plane Quantities
  H_hxfp = new TH1F("H_hxfp", Form("%s  X_{fp}; X_{fp} [cm]; Counts ", h_arm_name.Data()), hxfp_nbins, hxfp_xmin, hxfp_xmax);
  H_hyfp = new TH1F("H_hyfp", Form("%s  Y_{fp}; Y_{fp} [cm]; Counts ", h_arm_name.Data()), hyfp_nbins, hyfp_xmin, hyfp_xmax);
  H_hxpfp = new TH1F("H_hxpfp", Form("%s  X'_{fp}; X'_{fp} [rad]; Counts ", h_arm_name.Data()), hxpfp_nbins, hxpfp_xmin, hxpfp_xmax );
  H_hypfp = new TH1F("H_hypfp", Form("%s  Y'_{fp}; Y'_{fp} [rad]; Counts ", h_arm_name.Data()), hypfp_nbins, hypfp_xmin, hypfp_xmax);

  //Hadron arm Reconstructed Quantities 
  H_hytar = new TH1F("H_hytar", Form("%s  Y_{tar}; Y_{tar} [cm]; Counts ", h_arm_name.Data()), hytar_nbins, hytar_xmin, hytar_xmax);
  H_hxptar = new TH1F("H_hxptar", Form("%s  X'_{tar}; X'_{tar} [rad]; Counts ", h_arm_name.Data()), hxptar_nbins, hxptar_xmin, hxptar_xmax);
  H_hyptar = new TH1F("H_hyptar", Form("%s  Y'_{tar}; Y'_{tar} [rad]; Counts ", h_arm_name.Data()), hyptar_nbins, hyptar_xmin, hyptar_xmax );
  H_hdelta = new TH1F("H_hdelta", Form("%s  Momentum Acceptance, #delta; #delta [%%]; Counts ", h_arm_name.Data()), hdelta_nbins, hdelta_xmin, hdelta_xmax);
  

  //Target Reconstruction (Hall Coord. System) 
  H_htar_x = new TH1F("H_htar_x", Form("Fast Raster (%s) x-Vertex (Lab) ; x-Vertex [cm]; Counts ", h_arm_name.Data()), tarx_nbins, tarx_xmin, tarx_xmax);
  H_htar_y = new TH1F("H_htar_y", Form("Fast Raster (%s) y_Vertex (Lab) ; y-Vertex [cm]; Counts ", h_arm_name.Data()), tary_nbins, tary_xmin, tary_xmax);
  H_htar_z = new TH1F("H_htar_z", Form("%s z-Vertex (Lab)               ; z-Vertex [cm]; Counts ", h_arm_name.Data()), tarz_nbins, tarz_xmin, tarz_xmax);
  H_etar_x = new TH1F("H_etar_x", Form("Fast Raster (%s) x-Vertex (Lab) ; x-Target [cm]; Counts ", e_arm_name.Data()), tarx_nbins, tarx_xmin, tarx_xmax);
  H_etar_y = new TH1F("H_etar_y", Form("Fast Raster (%s) y-Vertex (Lab) ; y-Target [cm]; Counts ", e_arm_name.Data()), tary_nbins, tary_xmin, tary_xmax);
  H_etar_z = new TH1F("H_etar_z", Form("%s z-Vertex (Lab)               ; z-Vertex [cm]; Counts ", e_arm_name.Data()), tarz_nbins, tarz_xmin, tarz_xmax);

  //difference in reaction vertex z (user-defined)
  H_ztar_diff = new TH1F("H_ztar_diff", "z-Vertex Difference; (HMS-SHMS) z-Vertex Difference [cm]; Counts ", ztar_diff_nbins, ztar_diff_xmin, ztar_diff_xmax);

  //HMS / SHMS Collimator
  H_hXColl = new TH1F("H_hXColl", Form("%s X Collimator; X-Collimator [cm]; Counts ", h_arm_name.Data()), hXColl_nbins, hXColl_xmin, hXColl_xmax);
  H_hYColl = new TH1F("H_hYColl", Form("%s Y Collimator; Y-Collimator [cm]; Counts ", h_arm_name.Data()), hYColl_nbins, hYColl_xmin, hYColl_xmax); 
  H_eXColl = new TH1F("H_eXColl", Form("%s X Collimator; X-Collimator [cm]; Counts ", e_arm_name.Data()), eXColl_nbins, eXColl_xmin, eXColl_xmax);                                                                             
  H_eYColl = new TH1F("H_eYColl", Form("%s Y Collimator; Y-Collimator [cm]; Counts ", e_arm_name.Data()), eYColl_nbins, eYColl_xmin, eYColl_xmax);        

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


  
  //---------------------------------------------------------------------
  //---------HISTOGRAM CATEGORY: RANDOM COIN. BACKGROUND-----------------
  //---------------------------------------------------------------------

  H_ep_ctime_rand   = new TH1F("H_ep_ctime_rand", "ep Coincidence Time; ep Coincidence Time [ns]; Counts ", coin_nbins, coin_xmin, coin_xmax);
  H_W_rand          = new TH1F("H_W_rand",        "Invariant Mass, W", W_nbins, W_xmin, W_xmax); 
  H_Q2_rand         = new TH1F("H_Q2_rand",       "4-Momentum Transfer, Q^{2}", Q2_nbins, Q2_xmin, Q2_xmax); 
  H_xbj_rand        = new TH1F("H_xbj_rand",      "x-Bjorken", X_nbins, X_xmin, X_xmax);  
  H_nu_rand         = new TH1F("H_nu_rand",       "Energy Transfer, #nu", nu_nbins, nu_xmin, nu_xmax); 
  H_q_rand          = new TH1F("H_q_rand",        "3-Momentum Transfer, |#vec{q}|", q_nbins, q_xmin, q_xmax);
  H_Em_rand         = new TH1F("H_Em_rand",     "Missing Energy", Em_nbins, Em_xmin, Em_xmax);   
  H_Em_nuc_rand     = new TH1F("H_Em_nuc_rand", "Nuclear Missing Energy", Em_nuc_nbins, Em_nuc_xmin, Em_nuc_xmax); 
  H_Pm_rand         = new TH1F("H_Pm_rand",     "Missing Momentum, P_{miss}", Pm_nbins, Pm_xmin, Pm_xmax);
  H_MM_rand         = new TH1F("H_MM_rand",     "Missing Mass, M_{miss}", MM_nbins, MM_xmin, MM_xmax);        
  H_thxq_rand       = new TH1F("H_thxq_rand",   "In-Plane Angle, #theta_{pq}", thxq_nbins, thxq_xmin, thxq_xmax);
  H_thrq_rand       = new TH1F("H_thrq_rand",   "In-Plane Angle, #theta_{rq}", thrq_nbins, thrq_xmin, thrq_xmax);
  
  rand_HList->Add( H_ep_ctime_rand );
  rand_HList->Add( H_W_rand        );
  rand_HList->Add( H_Q2_rand       );
  rand_HList->Add( H_xbj_rand      );
  rand_HList->Add( H_nu_rand       );
  rand_HList->Add( H_q_rand        );
  rand_HList->Add( H_Em_rand       );
  rand_HList->Add( H_Em_nuc_rand   );
  rand_HList->Add( H_Pm_rand       );
  rand_HList->Add( H_MM_rand       );
  rand_HList->Add( H_thxq_rand     );
  rand_HList->Add( H_thrq_rand     );

  //--------------------------------------------------------------------------------
  //---------HISTOGRAM CATEGORY: RANDOM-SUBTRACTED COIN. BACKGROUND-----------------
  //--------------------------------------------------------------------------------

  H_ep_ctime_rand_sub   = new TH1F("H_ep_ctime_rand_sub", "ep Coincidence Time; ep Coincidence Time [ns]; Counts ", coin_nbins, coin_xmin, coin_xmax);
  H_W_rand_sub          = new TH1F("H_W_rand_sub",        "Invariant Mass, W", W_nbins, W_xmin, W_xmax); 
  H_Q2_rand_sub         = new TH1F("H_Q2_rand_sub",       "4-Momentum Transfer, Q^{2}", Q2_nbins, Q2_xmin, Q2_xmax); 
  H_xbj_rand_sub        = new TH1F("H_xbj_rand_sub",      "x-Bjorken", X_nbins, X_xmin, X_xmax);  
  H_nu_rand_sub         = new TH1F("H_nu_rand_sub",       "Energy Transfer, #nu", nu_nbins, nu_xmin, nu_xmax); 
  H_q_rand_sub          = new TH1F("H_q_rand_sub",        "3-Momentum Transfer, |#vec{q}|", q_nbins, q_xmin, q_xmax);
  H_Em_rand_sub         = new TH1F("H_Em_rand_sub",     "Missing Energy", Em_nbins, Em_xmin, Em_xmax);   
  H_Em_nuc_rand_sub     = new TH1F("H_Em_nuc_rand_sub", "Nuclear Missing Energy", Em_nuc_nbins, Em_nuc_xmin, Em_nuc_xmax); 
  H_Pm_rand_sub         = new TH1F("H_Pm_rand_sub",     "Missing Momentum, P_{miss}", Pm_nbins, Pm_xmin, Pm_xmax);
  H_MM_rand_sub         = new TH1F("H_MM_rand_sub",     "Missing Mass, M_{miss}", MM_nbins, MM_xmin, MM_xmax);        
  H_thxq_rand_sub       = new TH1F("H_thxq_rand_sub",   "In-Plane Angle, #theta_{pq}", thxq_nbins, thxq_xmin, thxq_xmax);
  H_thrq_rand_sub       = new TH1F("H_thrq_rand_sub",   "In-Plane Angle, #theta_{rq}", thrq_nbins, thrq_xmin, thrq_xmax);
  
  randSub_HList->Add( H_ep_ctime_rand_sub );
  randSub_HList->Add( H_W_rand_sub        );
  randSub_HList->Add( H_Q2_rand_sub       );
  randSub_HList->Add( H_xbj_rand_sub      );
  randSub_HList->Add( H_nu_rand_sub       );
  randSub_HList->Add( H_q_rand_sub        );
  randSub_HList->Add( H_Em_rand_sub       );
  randSub_HList->Add( H_Em_nuc_rand_sub   );
  randSub_HList->Add( H_Pm_rand_sub       );
  randSub_HList->Add( H_MM_rand_sub       );
  randSub_HList->Add( H_thxq_rand_sub     );
  randSub_HList->Add( H_thrq_rand_sub     );

  
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

  //Calculate Scaler Trigger Rates (EDTM subtracted already)
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

  
  if(analyze_data==true)
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

  else if(analyze_data==false)
    {
      cout << "SIMC ANALYSIS C++ CODE HAS NOT BEEN DONE YET ! ! !" << endl;
    } //END SIMC SET BRANCH ADDRESS


}

//_______________________________________________________________________________
Double_t baseAnalyzer::GetCoinTimePeak()
{
 
  cout << "Calling GetCoinTimePeak() . . . " <<  endl;
  
  // coin. time offset param (i.e., coin time peak value)
  Double_t ctime_offset = 0.0;
  
  // declare histogram to fill sample coin. time 
  TH1F *ctime_peak = new TH1F("ctime_peak", "Coin. Time Peak ", 200,-100,100);
  
  for(int ientry=0; ientry<10000; ientry++)
    {	  
      tree->GetEntry(ientry);
      // Fill sample histo to find peak
      ctime_peak->Fill(epCoinTime);	  	  

      cout << "SampleEventLoop: " << std::setprecision(2) << double(ientry) / 10000. * 100. << "  % " << std::flush << "\r";

    }
  
  
  // bin number corresponding to maximum bin content
  int binmax = ctime_peak->GetMaximumBin();
  
  // x-value corresponding to bin number with max content (i.e., peak)
  double xmax = ctime_peak->GetXaxis()->GetBinCenter(binmax);
  ctime_offset =  xmax;      

  cout << "coin time offset [ns] = " << ctime_offset << endl;  
  return ctime_offset; // in ns
  
   
}

//_______________________________________________________________________________
void baseAnalyzer::EventLoop()
{
  gROOT->SetBatch(1);
  cout << "Calling Base EventLoop() . . . " << endl;

  //Call Method to Set Collimator Graphical Cuts (In case it is used)
  CollimatorStudy();
  
  //Loop over Events
  
  if(analyze_data==true)
    {

      // Get Coin. Time peak to apply as an offset to center the coin. time peak at 0 ns    
      Double_t ctime_offset = GetCoinTimePeak();
	
      cout << "Loop over Data Events | nentries -->  " << nentries << endl;

      for(int ientry=0; ientry<nentries; ientry++)
	{
	  
	  tree->GetEntry(ientry);

	  
  
	  //cout << "ientry = " << ientry << endl;

	  //--------------CALCULATED KINEMATICS VARIABLES (IF THEY ARE NOT ALREADY DONE IN HCANA)-----------

	  th_x = xangle - th_e;  //hadron arm central angle for each particle
	  MM2 = MM*MM;           //Missing Mass Squared
 	  ztar_diff = htar_z - etar_z;  //reaction vertex z difference
	  
	  
	  // Calculate special missing energy to cut on background @ SRC kinematics (only for online analysis) Em = nu - Tp - T_n (for A>2 nuclei)	 
	  Em_src = nu - Tx - (sqrt(MN*MN + Pm*Pm) - MN); // assume kinetic energy of recoil system is that of a spectator SRC nucleon 
	  //cout << "Em_src = " << Em_src << endl;
	  //cout << "nu = " << nu << endl;
	  //cout << "Tx = " << Tx << endl;
	  //cout << "Pm = " << Pm << endl;
	  //cout << "MN = " << MN << endl;
	  
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

	  // CUTS (SPECIFIC TO DATA)
	 
	  // -- proton coincidence time cut ----
	  if(ePctime_cut_flag) {

	    // main coincidence time window cut
	    eP_ctime_cut = (epCoinTime-ctime_offset) >= ePctime_cut_min && (epCoinTime-ctime_offset) <= ePctime_cut_max;

	    // accidental coincidence (left/right of main coin. peak selected) as samples
	    eP_ctime_cut_rand_L = (epCoinTime-ctime_offset) >= ePctime_cut_max_L && (epCoinTime-ctime_offset) <= ePctime_cut_min_L ;
	    eP_ctime_cut_rand_R = (epCoinTime-ctime_offset) >= ePctime_cut_min_R && (epCoinTime-ctime_offset) <= ePctime_cut_max_R ;	    
	    
	    eP_ctime_cut_rand =  eP_ctime_cut_rand_L || eP_ctime_cut_rand_R;


	  }
	  else{
	    eP_ctime_cut=1;
	    eP_ctime_cut_rand =1;
	  }
	  
	  //----PID Cuts---- 
	  //SHMS calorimeter total normalized track energy
	  if(petot_trkNorm_pidCut_flag) {cpid_petot_trkNorm = pcal_etottracknorm>=cpid_petot_trkNorm_min && pcal_etottracknorm<=cpid_petot_trkNorm_max;}
	  else{cpid_petot_trkNorm=1;}

	  //SHMS Noble Gas Cherenkov
	  if(pngcer_pidCut_flag) {cpid_pngcer_NPE_Sum = pngcer_npesum>=cpid_pngcer_npeSum_min && pngcer_npesum<=cpid_pngcer_npeSum_max;}
	  else{cpid_pngcer_NPE_Sum=1;}

	  //SHMS Heavy Gas Cherenkov
	  if(phgcer_pidCut_flag) {cpid_phgcer_NPE_Sum = phgcer_npesum>=cpid_phgcer_npeSum_min && phgcer_npesum<=cpid_phgcer_npeSum_max;}
	  else{cpid_phgcer_NPE_Sum=1;}

	  c_pidCuts_shms = cpid_petot_trkNorm && cpid_pngcer_NPE_Sum && cpid_phgcer_NPE_Sum;
	  
	  //HMS calorimeter total normalized track energy
	  if(hetot_trkNorm_pidCut_flag) {cpid_hetot_trkNorm = hcal_etottracknorm>=cpid_hetot_trkNorm_min && hcal_etottracknorm<=cpid_hetot_trkNorm_max;}
	  else{cpid_hetot_trkNorm=1;}

	  //HMS Gas Cherenkov
	  if(hcer_pidCut_flag) {cpid_hcer_NPE_Sum = hcer_npesum>=cpid_hcer_npeSum_min && hcer_npesum<=cpid_hcer_npeSum_max;}
	  else{cpid_hcer_NPE_Sum=1;}

	  c_pidCuts_hms = cpid_hetot_trkNorm && cpid_hcer_NPE_Sum;

	  // combined hms/shms pid cuts 
	  c_pidCuts = c_pidCuts_shms && c_pidCuts_hms;


	   //----Acceptance Cuts----

	  // hadron arm
	  if(hdelta_cut_flag){c_hdelta = h_delta>=c_hdelta_min && h_delta<=c_hdelta_max;} 
	  else{c_hdelta=1;}

	  if(hxptar_cut_flag){c_hxptar = h_xptar>=c_hxptar_min && h_xptar<=c_hxptar_max;} 
	  else{c_hxptar=1;}

	  if(hyptar_cut_flag){c_hyptar = h_yptar>=c_hyptar_min && h_yptar<=c_hyptar_max;} 
	  else{c_hyptar=1;}

	  //Collimator CUTS
	  if(hmsCollCut_flag)  { hmsColl_Cut =  hms_Coll_gCut->IsInside(hYColl, hXColl);}
	  else{hmsColl_Cut=1;}
	  
	  c_accpCuts_hms = c_hdelta && c_hxptar && c_hyptar && hmsColl_Cut;
	  
	  // electron arm
	  if(edelta_cut_flag){c_edelta = e_delta>=c_edelta_min && e_delta<=c_edelta_max;} 
	  else{c_edelta=1;} 

	  if(exptar_cut_flag){c_exptar = e_xptar>=c_exptar_min && e_xptar<=c_exptar_max;} 
	  else{c_exptar=1;}

	  if(eyptar_cut_flag){c_eyptar = e_yptar>=c_eyptar_min && e_yptar<=c_eyptar_max;} 
	  else{c_eyptar=1;}

	  //Collimator Cuts
	  if(shmsCollCut_flag) { shmsColl_Cut =  shms_Coll_gCut->IsInside(eYColl, eXColl);}
	  else{shmsColl_Cut=1;}
	  
	  c_accpCuts_shms = c_edelta && c_exptar && c_eyptar && shmsColl_Cut;
  
	  // z-reaction vertex difference
	  if(ztarDiff_cut_flag){c_ztarDiff = ztar_diff>=c_ztarDiff_min && ztar_diff<=c_ztarDiff_max;} 
	  else{c_ztarDiff=1;}

	  // combined hms/shms acceptance cuts 
	  c_accpCuts = c_accpCuts_hms && c_accpCuts_shms && c_ztarDiff;
	  
	  //----Specialized Kinematics Cuts----

	  // H(e,e'p) Kinematics
	  
	  //Q2
	  if(Q2_heep_cut_flag){c_heep_Q2 = Q2>=c_heep_Q2_min && Q2<=c_heep_Q2_max;}
	  else{c_heep_Q2=1;}

	  //xbj
	  if(xbj_heep_cut_flag){c_heep_xbj = X>=c_heep_xbj_min && X<=c_heep_xbj_max;}
	  else{c_heep_xbj=1;}
	  
	  //Missing Energy, Em
	  if(Em_heep_cut_flag){c_heep_Em = Em>=c_heep_Em_min && Em<=c_heep_Em_max;}
	  else{c_heep_Em=1;}

	  //Invariant Mass, W
	  if(W_heep_cut_flag){c_heep_W = W>=c_heep_W_min && W<=c_heep_W_max;}
	  else{c_heep_W=1;}

	  //Missing Mass, MM = sqrt( E_recoil^2 - P_miss ^2 )
	  if(MM_heep_cut_flag){c_heep_MM = MM>=c_heep_MM_min && MM<=c_heep_MM_max;}
	  else{c_heep_MM=1;}


	  // H(e,e'p) singles ( e- trigger only)
	  c_kinHeepSing_Cuts = c_heep_Q2 && c_heep_W && c_heep_xbj;

	  // H(e,e'p) coin ( e- + p coin. trigger )
	  c_kinHeepCoin_Cuts = c_heep_Q2 && c_heep_xbj && c_heep_Em && c_heep_W && c_heep_MM;
	  
	  // CaFe A(e,e'p) Mean-Field (MF) Kinematic Cuts

	  // Q2
	  if(Q2_MF_cut_flag){c_MF_Q2 = Q2>=c_MF_Q2_min && Q2<=c_MF_Q2_max;}
	  else{c_MF_Q2=1;}

	  // Pm
	  if(Pm_MF_cut_flag){c_MF_Pm = Pm>=c_MF_Pm_min && Pm<=c_MF_Pm_max;}
	  else{c_MF_Pm=1;}

	  // Em ( require this cut ONLY for deuteron)
	  if(Em_d2MF_cut_flag && tgt_type=="LD2"){c_d2MF_Em = Em_nuc>=c_d2MF_Em_min && Em_nuc <= c_d2MF_Em_max;}
	  else{c_d2MF_Em=1;}

	  // Em ( require this cut ONLY for A>2 nuclei)
	  if(Em_MF_cut_flag && tgt_type!="LD2"){c_MF_Em = Em_nuc >= c_MF_Em_min && Em_nuc <= c_MF_Em_max;}
	  else{c_MF_Em=1;}
	  
	  c_kinMF_Cuts = c_MF_Q2 && c_MF_Pm && c_d2MF_Em && c_MF_Em;
	    
	  // CaFe A(e,e'p) Short-Range Correlations (SRC) Kinematic Cuts

	  // Q2
	  if(Q2_SRC_cut_flag){c_SRC_Q2 = Q2>=c_SRC_Q2_min && Q2<=c_SRC_Q2_max;}
	  else{c_SRC_Q2=1;}

	  // Pm
	  if(Pm_SRC_cut_flag){c_SRC_Pm = Pm>=c_SRC_Pm_min && Pm<=c_SRC_Pm_max;}
	  else{c_SRC_Pm=1;}

	  // Xbj
	  if(Xbj_SRC_cut_flag){c_SRC_Xbj = X >= c_SRC_Xbj_min && X <= c_SRC_Xbj_max;}
	  else{c_SRC_Xbj=1;}
	  
	  // theta_rq
	  if(thrq_SRC_cut_flag){c_SRC_thrq = th_rq >= c_SRC_thrq_min && th_rq <= c_SRC_thrq_max;}
	  else{c_SRC_thrq=1;}
	  
	  // Em ( require this cut ONLY for deuteron)
	  if(Em_d2SRC_cut_flag && tgt_type=="LD2"){c_d2SRC_Em = Em_nuc>=c_d2SRC_Em_min && Em_nuc <= c_d2SRC_Em_max;}
	  else{c_d2SRC_Em=1;}

	  // Em ( require this cut ONLY for A>2 nuclei)
	  if(Em_SRC_cut_flag && tgt_type!="LD2"){c_SRC_Em = Em_src>0. && Em_nuc <= Em_src; // put lower bound on Em_src cut
	    //cout << "c_SRC_Em = " << c_SRC_Em << endl;
	  }
	  else{c_SRC_Em=1;}
	  

	  c_kinSRC_Cuts = c_SRC_Q2 && c_SRC_Pm && c_SRC_Xbj && c_SRC_thrq && c_d2SRC_Em && c_SRC_Em;


	 
			  	 
	  // ----- Combine All CUTS -----

	  // user pre-determined analysis kinematics cuts

	  if(analysis_cut=="lumi"){ 
	    c_baseCuts =  e_delta>=-10. && e_delta<=22. && c_pidCuts_shms;
	  }
	  else if(analysis_cut=="optics"){  // will need to call Holly's script that generates optics plots (from raw ROOTfile)
	    c_baseCuts =  c_pidCuts_shms;
	  }
	  else if(analysis_cut=="heep_singles"){
	    c_baseCuts =  c_accpCuts_shms && c_pidCuts_shms && c_kinHeepSing_Cuts && (eP_ctime_cut=1); //setting eP_ctime=1 guarantees cut will always pass (i.e. turned coin time cut OFF)
	  }
	  else if(analysis_cut=="heep_coin"){
	    c_baseCuts =  c_accpCuts && c_pidCuts && c_kinHeepCoin_Cuts;
	  }
	  else if(analysis_cut=="MF"){
	    c_baseCuts =  c_accpCuts && c_pidCuts && c_kinMF_Cuts;
	  }
	  else if(analysis_cut=="SRC"){
	    c_baseCuts =  c_accpCuts && c_pidCuts && c_kinSRC_Cuts;
	  }
	  
	  
	  
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
	      
	      //cout << "passed BCM Cut !" << endl;
	      bool event_type_cut = false;
	      
	      if( (analysis_cut=="heep_singles") || (analysis_cut=="lumi") || (analysis_cut=="optics") || (analysis_cut=="bcm_calib") ){
		event_type_cut = (gevtyp==1 || gevtyp==3);  // use this to calculate live time for shms singles events only                               
              }
	      else if((analysis_cut=="heep_coin") || (analysis_cut=="MF") || (analysis_cut=="SRC")){
		event_type_cut = (gevtyp == 4);} //use this to calculate live time for coin. events only                                                                          
	      //Count Accepted EDTM events (With bcm current cut: to be used in total edtm live time calculation)
	      if(c_edtm && event_type_cut){ total_edtm_accp_bcm_cut++;}
	      
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
		  //cout << "passed NO EDTM Cut !" << endl;   
		  //Calculate HMS Tracking Efficiency Components
		  if(good_hms_did){ h_did++;}
		  if(good_hms_should){ h_should++; }
		  
		  //Calculate SHMS Tracking Efficiency Components
		  if(good_shms_did){ p_did++;}
		  if(good_shms_should){ p_should++; }

		  
		  //----------------------Fill DATA Histograms-----------------------

		  
		  //2D Kin plots to help clean out online Em data
		  if(c_accpCuts && c_pidCuts && eP_ctime_cut){
		    H_Em_nuc_vs_Pm ->Fill(Pm, Em_nuc);
		    H_Em_src_vs_Pm ->Fill(Pm, Em_src);
		  }
		  
		  if(c_baseCuts){
		    
		    
		    // full coin. time spectrum with all other cuts  
		    H_ep_ctime_total->Fill(epCoinTime-ctime_offset); 
		  
		    
		    // select "TRUE COINCIDENCE " (electron-proton from same "beam bunch" form a coincidence)
		    if(eP_ctime_cut)
		      {
			
			//--------------------------------------------------------------------
			//---------HISTOGRAM CATEGORY: Particle Identification (PID)----------
			//--------------------------------------------------------------------
			
			//Coincidence Time		      
			H_ep_ctime->Fill(epCoinTime-ctime_offset); // fill coin. time and apply the offset
			
			//Fill HMS Detectors
			H_hCerNpeSum->Fill(hcer_npesum);
			H_hCalEtotNorm->Fill(hcal_etotnorm);
			H_hCalEtotTrkNorm->Fill(hcal_etottracknorm);
			H_hHodBetaNtrk->Fill(hhod_beta_ntrk);
			H_hHodBetaTrk->Fill(hhod_gtr_beta);
			
			//Fill SHMS Detectors
			H_pNGCerNpeSum->Fill(pngcer_npesum);
			H_pHGCerNpeSum->Fill(phgcer_npesum);
			H_pCalEtotNorm->Fill(pcal_etotnorm);
			H_pCalEtotTrkNorm->Fill(pcal_etottracknorm);
			H_pHodBetaNtrk->Fill(phod_beta_ntrk);
			H_pHodBetaTrk->Fill(phod_gtr_beta);
			
			//Fill 2D PID Correlations
			H_hcal_vs_hcer->Fill(hcal_etottracknorm, hcer_npesum);
			H_pcal_vs_phgcer->Fill(pcal_etottracknorm, phgcer_npesum);  
			H_pcal_vs_pngcer->Fill(pcal_etottracknorm, pngcer_npesum);   
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
			H_Em_src   ->Fill(Em_src);
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
			
			
		      } // end TRUE COINCIDENCE time cut
		    
		    
		    // select "ACCIDENTAL COINCIDENCE BACKGROUND" left/right of main coin. peak as a sample to estimate background underneath main coin. peak
		    // background underneath main peak: electron-proton from same "beam bunch" form a random ("un-correlated") coincidence

		    if(ePctime_cut_flag && eP_ctime_cut_rand)
		      
		      {
			// Only histograms of selected variables of interest will be filled with background
			
			H_ep_ctime_rand->  Fill ( epCoinTime-ctime_offset );
			H_W_rand       ->  Fill (W);       
			H_Q2_rand      ->  Fill (Q2);      
			H_xbj_rand     ->  Fill (X);     
			H_nu_rand      ->  Fill (nu);      
			H_q_rand       ->  Fill (q);       
			H_Em_rand      ->  Fill (Em);      
			H_Em_nuc_rand  ->  Fill (Em_nuc);  
			H_Pm_rand      ->  Fill (Pm);      
			H_MM_rand      ->  Fill (MM);      
			H_thxq_rand    ->  Fill (th_xq/dtr);    
			H_thrq_rand    ->  Fill (ph_xq/dtr);    
				       	
		      }  		    		   
		    
		   
		    
		  }  //----------------------END: Fill DATA Histograms-----------------------		  		 		  
		  
		  
		} //------END: REQUIRE "NO EDTM" CUT TO FILL DATA HISTOGRAMS-----
	      
	    }  //-----END: BCM Current Cut------
	  
	  
	  
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

  if(analyze_data==false)
    {
      cout << "SIMC ANALYSIS needs to be done . . . " << endl;
    }
  
}

//_______________________________________________________________________________
void baseAnalyzer::RandSub()
{
  cout << "Calling RandSub() " << endl;

  /*
    Brief: This methods carries out the subtraction of random coincidences (outside coin peak selection) 
    from real coincidences (within coin peak selected) for various histograms 
  */

  // Scale Down (If necessary) the randoms before subtracting it from the reals
  // NOTE: the accidentals selected to the left/right of main coin. peak have a window width of dt_acc_L and dt_acc_R
  // and therefore the relevant histograms with accidentals selection  must be properly scaled down to the main coin. peak window
  // as follows:  scale_factor = dt_coin_peak / ( dt_acc_L + dt_acc_R ) --> ratio of main coin. time window width to (sum of accidental window width left/right of main peak)
  // the random coincidences can be scaled down by the following factor:

  // main. coin. peak is assumed to have been centered at zero
  dt_coin_peak =  ePctime_cut_max - ePctime_cut_min;
  dt_acc_R = abs(ePctime_cut_max_R - ePctime_cut_min_R);
  dt_acc_L = abs(ePctime_cut_max_L - ePctime_cut_min_L);
  
  P_scale_factor = dt_coin_peak /  (dt_acc_L + dt_acc_R);

  cout << "coin_peak_window_width [ns] = "         << dt_coin_peak << endl;
  cout << "accidental_width_LEFT [ns] = " << dt_acc_L << endl;
  cout << "accidental_width_RIGHT [ns] = " << dt_acc_R << endl;  
  cout << "P_scale_factor = "         << P_scale_factor << endl;


  //----Scale Down the random coincidences histograms-----
  // ----(other than the coin. histograms themselves)----

  H_W_rand       ->  Scale( P_scale_factor );
  H_Q2_rand      ->  Scale( P_scale_factor );
  H_xbj_rand     ->  Scale( P_scale_factor );
  H_nu_rand      ->  Scale( P_scale_factor );
  H_q_rand       ->  Scale( P_scale_factor );
  H_Em_rand      ->  Scale( P_scale_factor );
  H_Em_nuc_rand  ->  Scale( P_scale_factor );
  H_Pm_rand      ->  Scale( P_scale_factor );
  H_MM_rand      ->  Scale( P_scale_factor );
  H_thxq_rand    ->  Scale( P_scale_factor );
  H_thrq_rand    ->  Scale( P_scale_factor );

  
  // -----Carry out the randoms subtraction------
  H_W_rand_sub       -> Add(H_W      ,H_W_rand      , 1, -1);
  H_Q2_rand_sub      -> Add(H_Q2     ,H_Q2_rand     , 1, -1);
  H_xbj_rand_sub     -> Add(H_xbj    ,H_xbj_rand    , 1, -1);
  H_nu_rand_sub      -> Add(H_nu     ,H_nu_rand     , 1, -1);
  H_q_rand_sub       -> Add(H_q      ,H_q_rand      , 1, -1);
  H_Em_rand_sub      -> Add(H_Em     ,H_Em_rand     , 1, -1);
  H_Em_nuc_rand_sub  -> Add(H_Em_nuc ,H_Em_nuc_rand , 1, -1);
  H_Pm_rand_sub      -> Add(H_Pm     ,H_Pm_rand     , 1, -1);
  H_MM_rand_sub      -> Add(H_MM     ,H_MM_rand     , 1, -1);
  H_thxq_rand_sub    -> Add(H_thxq   ,H_thxq_rand   , 1, -1);
  H_thrq_rand_sub    -> Add(H_thrq   ,H_thrq_rand   , 1, -1);  

      
  // Get Counts of "good events for saving to CaFe Report File"
  total_bins = H_W->GetNbinsX();  //Get total number of bins (excluding overflow) (same for total, reals randoms, provied same histo range)
  
  W_total = H_W          ->IntegralAndError(1, total_bins, W_total_err);
  W_real  = H_W_rand_sub ->IntegralAndError(1, total_bins, W_real_err);
  W_rand  = H_W_rand     ->IntegralAndError(1, total_bins, W_rand_err);
  cout << Form("W_total = %.3f", W_total) << endl;
  cout << Form("W_real = %.3f", W_real) << endl;
  cout << Form("W_rand = %.3f", W_rand) << endl;

  W_total_rate = W_total / total_time_bcm_cut;  // # good elastic proton event rate
  W_real_rate = W_real / total_time_bcm_cut;
   
  total_bins = H_Pm->GetNbinsX(); 
  Pm_total = H_Pm          ->IntegralAndError(1, total_bins, Pm_total_err);
  Pm_real  = H_Pm_rand_sub ->IntegralAndError(1, total_bins, Pm_real_err);
  Pm_rand  = H_Pm_rand     ->IntegralAndError(1, total_bins, Pm_rand_err);

  Pm_real_rate = Pm_real / total_time_bcm_cut;

  
  total_bins = H_Em->GetNbinsX(); 
  Em_total = H_Em          ->IntegralAndError(1, total_bins, Em_total_err);
  Em_real  = H_Em_rand_sub ->IntegralAndError(1, total_bins, Em_real_err);
  Em_rand  = H_Em_rand     ->IntegralAndError(1, total_bins, Em_rand_err);
  
  total_bins = H_Em_nuc->GetNbinsX(); 
  Em_nuc_total = H_Em_nuc          ->IntegralAndError(1, total_bins, Em_nuc_total_err);
  Em_nuc_real  = H_Em_nuc_rand_sub ->IntegralAndError(1, total_bins, Em_nuc_real_err);
  Em_nuc_rand  = H_Em_nuc_rand     ->IntegralAndError(1, total_bins, Em_nuc_rand_err);
  
  total_bins = H_MM->GetNbinsX(); 
  MM_total = H_MM          ->IntegralAndError(1, total_bins, MM_total_err);
  MM_real  = H_MM_rand_sub ->IntegralAndError(1, total_bins, MM_real_err);
  MM_rand  = H_MM_rand     ->IntegralAndError(1, total_bins, MM_rand_err);
  
  
  
}

//_________________________________________________________
void baseAnalyzer::CollimatorStudy()
{

  //Method to study various collimator cuts on the H(e,e'p) and D(e,e'p)n  Yield across Ytar, Y'tar, X'tar and delta

  cout << "Calling CollimatorStudy() . . . " << endl;
  
  //Scaling the HMS/SHMS Collimator Cuts
  hms_hsize = hms_scale*hms_hsize;  //The scale factor is read from set_heep_cuts.inp
  hms_vsize = hms_scale*hms_vsize;
  
  shms_hsize = shms_scale*shms_hsize;
  shms_vsize = shms_scale*shms_vsize;  

  //Define HMS Collimator Shape
  hms_Coll_gCut = new TCutG("hmsCollCut", 8 );
  hms_Coll_gCut->SetVarX("X");
  hms_Coll_gCut->SetVarY("Y");
 
  hms_Coll_gCut->SetPoint(0,  hms_hsize,     hms_vsize/2.);
  hms_Coll_gCut->SetPoint(1,  hms_hsize/2.,  hms_vsize   );
  hms_Coll_gCut->SetPoint(2, -hms_hsize/2.,  hms_vsize   );
  hms_Coll_gCut->SetPoint(3, -hms_hsize,     hms_vsize/2.);
  hms_Coll_gCut->SetPoint(4, -hms_hsize,    -hms_vsize/2.);
  hms_Coll_gCut->SetPoint(5, -hms_hsize/2., -hms_vsize   );
  hms_Coll_gCut->SetPoint(6,  hms_hsize/2., -hms_vsize   );
  hms_Coll_gCut->SetPoint(7,  hms_hsize,    -hms_vsize/2.);
  hms_Coll_gCut->SetPoint(8,  hms_hsize,     hms_vsize/2.);

  //Define SHMS Collimator Shape
  shms_Coll_gCut = new TCutG("shmsCollCut", 8 );
  shms_Coll_gCut->SetVarX("X");
  shms_Coll_gCut->SetVarY("Y");
 
  shms_Coll_gCut->SetPoint(0,  shms_hsize,     shms_vsize/2.);
  shms_Coll_gCut->SetPoint(1,  shms_hsize/2.,  shms_vsize   );
  shms_Coll_gCut->SetPoint(2, -shms_hsize/2.,  shms_vsize   );
  shms_Coll_gCut->SetPoint(3, -shms_hsize,     shms_vsize/2.);
  shms_Coll_gCut->SetPoint(4, -shms_hsize,    -shms_vsize/2.);
  shms_Coll_gCut->SetPoint(5, -shms_hsize/2., -shms_vsize   );
  shms_Coll_gCut->SetPoint(6,  shms_hsize/2., -shms_vsize   );
  shms_Coll_gCut->SetPoint(7,  shms_hsize,    -shms_vsize/2.);
  shms_Coll_gCut->SetPoint(8,  shms_hsize,     shms_vsize/2.);

  cout << "Ending CollimatorStudy() . . . " << endl;


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

  //Convert Scaler Trigger/EDTM Rates from Hz to kHz 
  S1XscalerRate_bcm_cut   = S1XscalerRate_bcm_cut   / 1000.;
  TRIG1scalerRate_bcm_cut = TRIG1scalerRate_bcm_cut / 1000.;
  TRIG2scalerRate_bcm_cut = TRIG2scalerRate_bcm_cut / 1000.;
  TRIG3scalerRate_bcm_cut = TRIG3scalerRate_bcm_cut / 1000.;
  TRIG4scalerRate_bcm_cut = TRIG4scalerRate_bcm_cut / 1000.;
  TRIG5scalerRate_bcm_cut = TRIG5scalerRate_bcm_cut / 1000.;
  TRIG6scalerRate_bcm_cut = TRIG6scalerRate_bcm_cut / 1000.;
  EDTMscalerRate_bcm_cut  = EDTMscalerRate_bcm_cut  / 1000.;

  /*---- Apply Pre-scale factor to accepted triggers (to make comparisons with scaler triggers)
  total_trig1_accp_bcm_cut = total_trig1_accp_bcm_cut * Ps1_factor;
  total_trig2_accp_bcm_cut = total_trig2_accp_bcm_cut * Ps2_factor;  
  total_trig3_accp_bcm_cut = total_trig3_accp_bcm_cut * Ps3_factor;  
  total_trig4_accp_bcm_cut = total_trig4_accp_bcm_cut * Ps4_factor;  
  total_trig5_accp_bcm_cut = total_trig5_accp_bcm_cut * Ps5_factor;  
  total_trig6_accp_bcm_cut = total_trig6_accp_bcm_cut * Ps6_factor;  
  */
  //Calculate Accepted Trigger/EDTM Rates in kHz
  TRIG1accpRate_bcm_cut = (total_trig1_accp_bcm_cut / total_time_bcm_cut ) / 1000.;
  TRIG2accpRate_bcm_cut = (total_trig2_accp_bcm_cut / total_time_bcm_cut ) / 1000.;
  TRIG3accpRate_bcm_cut = (total_trig3_accp_bcm_cut / total_time_bcm_cut ) / 1000.;
  TRIG4accpRate_bcm_cut = (total_trig4_accp_bcm_cut / total_time_bcm_cut ) / 1000.;
  TRIG5accpRate_bcm_cut = (total_trig5_accp_bcm_cut / total_time_bcm_cut ) / 1000.;
  TRIG6accpRate_bcm_cut = (total_trig6_accp_bcm_cut / total_time_bcm_cut ) / 1000.;
  EDTMaccpRate_bcm_cut  = (total_edtm_accp_bcm_cut  / total_time_bcm_cut ) / 1000.;
  
  //Calculate Pure Computer Live Time (numerator->accepted tdc trig requires NO EDTM :: denominator -> EDTM has already been subtracted from scaler counts)
  //Pre-Scale factor has been accounted 
  cpuLT_trig1 = total_trig1_accp_bcm_cut * Ps1_factor / (total_trig1_scaler_bcm_cut);
  cpuLT_trig2 = total_trig2_accp_bcm_cut * Ps2_factor / (total_trig2_scaler_bcm_cut);
  cpuLT_trig3 = total_trig3_accp_bcm_cut * Ps3_factor / (total_trig3_scaler_bcm_cut);
  cpuLT_trig4 = total_trig4_accp_bcm_cut * Ps4_factor / (total_trig4_scaler_bcm_cut);
  cpuLT_trig5 = total_trig5_accp_bcm_cut * Ps5_factor / (total_trig5_scaler_bcm_cut);
  cpuLT_trig6 = total_trig6_accp_bcm_cut * Ps6_factor / (total_trig6_scaler_bcm_cut);

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

  /*
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
  */
  
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

  //For testing purposes and online experiment production (do not scale by charge or det. inefficieny)
  FullWeight = 1.; // / total_charge_bcm_cut;

  if((analysis_cut=="heep_singles") || (analysis_cut=="optics") || (analysis_cut=="lumi") || (analysis_cut=="bcm_calib")){ // For CaFe, PS2_factor is pre-scale factor for SHMS EL-REAL
    FullWeight = Ps2_factor; // if accepted trigger pre-scaled, scale by pre-scale factor to recover events
  }
  else{ // else use pre-scale factor determined from trig_type input parameter
    FullWeight = Ps_factor; 
  }
  //Scale Data Histograms by Full Weight (Each run for a particular kinematics can then be combined, once they are scaled by the FullWeight)
  
  //----SCALE HISTOGRAMS BY LOOPING OVER LISTS----
  
  //determine what class types are in the list
  TString class_name;
  
  
  //-----------------------------------------------------
  //Lopp over pid_HList of histogram objects 
  //----------------------------------------------------
  for(int i=0; i<pid_HList->GetEntries(); i++) {
    //Get the class name for each element on the list (either "TH1F" or TH2F")
    class_name = pid_HList->At(i)->ClassName();
    //Read ith histograms in the list from current run
    if(class_name=="TH1F") {
      //Get and scale histogram from the list
      h_i = (TH1F *)pid_HList->At(i); h_i->Scale(FullWeight); 
    }
    if(class_name=="TH2F") {
      //Get and scale histogram from the list
      h2_i = (TH2F *)pid_HList->At(i); h2_i->Scale(FullWeight);
    }   
  }//end loop over pid_HList
  
  //-----------------------------------------------------
  //Lopp over kin_HList of histogram objects 
  //----------------------------------------------------
  for(int i=0; i<kin_HList->GetEntries(); i++) {
    //Get the class name for each element on the list (either "TH1F" or TH2F")
    class_name = kin_HList->At(i)->ClassName();
    //Read ith histograms in the list from current run
    if(class_name=="TH1F") {
      //Get and scale histogram from the list
      h_i = (TH1F *)kin_HList->At(i); h_i->Scale(FullWeight); 
    }
    if(class_name=="TH2F") {
      //Get and scale histogram from the list
      h2_i = (TH2F *)kin_HList->At(i); h2_i->Scale(FullWeight);
    }   
  }//end loop over kin_HList	
  
  //-----------------------------------------------------
  //Lopp over accp_HList of histogram objects 
  //----------------------------------------------------
  for(int i=0; i<accp_HList->GetEntries(); i++) {
    //Get the class name for each element on the list (either "TH1F" or TH2F")
    class_name = accp_HList->At(i)->ClassName();
    //Read ith histograms in the list from current run
    if(class_name=="TH1F") {
      //Get and scale histogram from the list
      h_i = (TH1F *)accp_HList->At(i); h_i->Scale(FullWeight); 
    }
    if(class_name=="TH2F") {
      //Get and scale histogram from the list
      h2_i = (TH2F *)accp_HList->At(i); h2_i->Scale(FullWeight);
    }   
  }//end loop over accp_HList


  //-----------------------------------------------------
  //Lopp over rand_HList of histogram objects 
  //----------------------------------------------------
  for(int i=0; i<rand_HList->GetEntries(); i++) {
    //Get the class name for each element on the list (either "TH1F" or TH2F")
    class_name = rand_HList->At(i)->ClassName();
    //Read ith histograms in the list from current run
    if(class_name=="TH1F") {
      //Get and scale histogram from the list
      h_i = (TH1F *)rand_HList->At(i); h_i->Scale(FullWeight); 
    }
    if(class_name=="TH2F") {
      //Get and scale histogram from the list
      h2_i = (TH2F *)rand_HList->At(i); h2_i->Scale(FullWeight);
    }   
  }//end loop over accp_HList
  

  //-----------------------------------------------------
  //Lopp over randSub_HList of histogram objects 
  //----------------------------------------------------
  for(int i=0; i<randSub_HList->GetEntries(); i++) {
    //Get the class name for each element on the list (either "TH1F" or TH2F")
    class_name = randSub_HList->At(i)->ClassName();
    //Read ith histograms in the list from current run
    if(class_name=="TH1F") {
      //Get and scale histogram from the list
      h_i = (TH1F *)randSub_HList->At(i); h_i->Scale(FullWeight); 
    }
    if(class_name=="TH2F") {
      //Get and scale histogram from the list
      h2_i = (TH2F *)randSub_HList->At(i); h2_i->Scale(FullWeight);
    }   
  }//end loop over accp_HList
  

  //Call the randoms subtraction method, provided there was a coin. time cut flag  (after scaling all histograms above)
  RandSub();
 
}

//_______________________________________________________________________________
void baseAnalyzer::WriteHist()
{
  /*
    Brief: Method to write histograms to a ROOTfile
  */

  cout << "Calling WriteHist() . . ." << endl;

  
  //Write Data Histograms
  if(analyze_data==true)
    {
      //Create Output ROOTfile
      outROOT = new TFile(data_OutputFileName, "RECREATE");

      //Make directories to store histograms based on category
      outROOT->mkdir("pid_plots");
      outROOT->mkdir("kin_plots");
      outROOT->mkdir("accp_plots");

      outROOT->mkdir("rand_plots");
      outROOT->mkdir("randSub_plots");

      //Write PID histos to pid_plots directory
      outROOT->cd("pid_plots");
      pid_HList->Write();
      
      //Write Kinematics histos to kin_plots directory
      outROOT->cd("kin_plots");
      kin_HList->Write();

      //Write Acceptance histos to accp_plots directory
      outROOT->cd("accp_plots");
      accp_HList->Write();
      
      //Write selected Random histos to rand_plots directory
      outROOT->cd("rand_plots");
      rand_HList->Write();
      
      //Write selected Random-Subtracted histos to randSub_plots directory
      outROOT->cd("randSub_plots");
      randSub_HList->Write();


      
      //Close File
      outROOT->Close();
    }

  
}
//_______________________________________________________________________________
void baseAnalyzer::WriteReport()
{
  
  /*  Method to write charge, efficiencies, live time and other relevant quantities to a data file
      on a run-by-run basis, meaning, each run that a report is written separately for each run.
   */
  
  cout << "Calling WriteReport() . . ." << endl;
  
  if(analyze_data==true){


    if( (analysis_cut=="MF") || (analysis_cut=="SRC")) {

      cafe_Ib_simc = stod(split(FindString("cafe_Ib_simc",    input_SIMCinfo_FileName.Data())[0], '=')[1]);

      total_simc_counts = stod(split(FindString(Form("%s_%s_counts", tgt_type.Data(), analysis_cut.Data()),    input_SIMCinfo_FileName.Data())[0], '=')[1]); // [counts]
      total_simc_time = stod(split(FindString(Form("%s_%s_time", tgt_type.Data(), analysis_cut.Data()),    input_SIMCinfo_FileName.Data())[0], '=')[1]); // [hr]
      simc_cafe_rates = total_simc_counts / (total_simc_time * 3600.); //[Hz]
      
      // [mC]                [uC / sec]        [hr]      [sec]/[hr]  0.001 mC / 1 uC
      total_simc_charge =  cafe_Ib_simc * total_simc_time * 3600. * 1e-3;  
    }
    
    else if( (analysis_cut=="heep_singles") || (analysis_cut=="heep_coin") ) {

      heep_Ib_simc = stod(split(FindString("heep_Ib_simc",    input_SIMCinfo_FileName.Data())[0], '=')[1]);

      heep_kin0_counts = stod(split(FindString("heep_kin0_counts",    input_SIMCinfo_FileName.Data())[0], '=')[1]); 
      heep_kin1_counts = stod(split(FindString("heep_kin1_counts",    input_SIMCinfo_FileName.Data())[0], '=')[1]); 
      heep_kin2_counts = stod(split(FindString("heep_kin2_counts",    input_SIMCinfo_FileName.Data())[0], '=')[1]); 

      heep_kin0_time = stod(split(FindString("heep_kin0_time",    input_SIMCinfo_FileName.Data())[0], '=')[1]); // [hr]
      heep_kin1_time = stod(split(FindString("heep_kin1_time",    input_SIMCinfo_FileName.Data())[0], '=')[1]); // [hr]
      heep_kin2_time = stod(split(FindString("heep_kin2_time",    input_SIMCinfo_FileName.Data())[0], '=')[1]); // [hr]

      heep_kin0_rates = heep_kin0_counts / (heep_kin0_time * 3600.); // [Hz]
      heep_kin1_rates = heep_kin1_counts / (heep_kin1_time * 3600.); // [Hz]
      heep_kin2_rates = heep_kin2_counts / (heep_kin2_time * 3600.); // [Hz]

      // [mC]                [uC / sec]        [hr]      [sec]/[hr]  0.001 mC / 1 uC
      heep_kin0_charge =  heep_Ib_simc * heep_kin0_time * 3600. * 1e-3;
      heep_kin1_charge =  heep_Ib_simc * heep_kin1_time * 3600. * 1e-3;
      heep_kin2_charge =  heep_Ib_simc * heep_kin2_time * 3600. * 1e-3;
    
    }
    
    
    
    //Check if file already exists
    in_file.open(output_ReportFileName.Data());

    if(!in_file.fail()){
      cout << Form("Report File for run %d exists, will overwrite it . . . ", run) << endl;
    }    
    else if(in_file.fail()){
      cout << "Report File does NOT exist, will create one . . . " << endl;
    }
    
    out_file.open(output_ReportFileName);
    out_file << Form("# Run %d Data Analysis Summary", run)<< endl;
    out_file << "                                     " << endl;
    out_file << "# =:=:=:=:=:=:=:=:=:=:=:=:=:=:" << endl;
    out_file << "# General Run Configuration                              " << endl;
    out_file << "# =:=:=:=:=:=:=:=:=:=:=:=:=:=:" << endl;
    out_file << "                                     " << endl;
    out_file << Form("run_number: %d                     ", run) << endl;
    out_file << "" << endl;
    out_file << Form("start_of_run = %s", start_of_run.Data())<< endl;
    out_file << Form("end_of_run = %s", end_of_run.Data())<< endl;
    out_file << "" << endl;    
    out_file << Form("kin_type: %s                     ", analysis_cut.Data()) << endl;
    out_file << Form("daq_mode: %s                     ", daq_mode.Data()) << endl;
    out_file << Form("events_replayed: %lld              ", nentries ) << endl;
    out_file << "" << endl;
    out_file << Form("beam_energy [GeV]: %.4f          ", beam_energy ) << endl;          
    out_file << Form("target_name: %s                       ", tgt_type.Data() ) << endl;
    out_file << Form("target_amu: %.6f                 ", tgt_mass        ) << endl;      
    out_file << "" << endl;      
    out_file << Form("hms_h_particle_mass [GeV]: %.6f          ",  hms_part_mass ) << endl;          
    out_file << Form("hms_h_momentum [GeV/c]: %.4f             ",  hms_p ) << endl;
    out_file << Form("hms_h_angle [deg]: %.4f                  ",  hms_angle ) << endl;          
    out_file << "" << endl;      
    out_file << Form("shms_e_particle_mass [GeV]: %.6f          ",  shms_part_mass ) << endl;          
    out_file << Form("shms_e_momentum [GeV/c]: %.4f             ",  shms_p ) << endl;
    out_file << Form("shms_e_angle [deg]: %.4f                  ",  shms_angle ) << endl;  
    out_file << "" << endl;      
    out_file << Form("%s_Current_Threshold [uA]: >%.2f ", bcm_type.Data(), bcm_thrs) << endl;
    out_file << Form("beam_on_target [sec]: %.3f       ", total_time_bcm_cut) << endl;
    out_file << Form("%s_Average_Current [uA]: %.3f ", bcm_type.Data(), avg_current_bcm_cut ) << endl;
    out_file << Form("%s_Charge [mC]: %.3f ", bcm_type.Data(), total_charge_bcm_cut ) << endl;
    out_file << "" << endl;
    if(analysis_cut=="heep_singles")
      {
	out_file << "# =:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:" << endl;
	out_file << "# CaFe H(e,e')p  Singles Counts    " << endl;
	out_file << "# =:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:" << endl;
    
	out_file << Form("heep_total_singles_counts  : %.3f ", W_total) << endl;
	out_file << Form("heep_total_singles_rate [Hz]  : %.3f ", W_total_rate) << endl;
	out_file << "" << endl;
	// check if shms angle is within 0.1 deg of nominal
	if(abs(8.3-shms_angle)<0.1) { out_file << Form("simc_heep_kin0_rates (shms=8.3 deg) [Hz] x (%.1f uA/%.1f uA) : %.3f ", avg_current_bcm_cut, heep_Ib_simc,  heep_kin0_rates * (avg_current_bcm_cut/heep_Ib_simc) ) << endl;}
	if(abs(7.5-shms_angle)<0.1) { out_file << Form("simc_heep_kin1_rates (shms=7.5 deg) [Hz] x (%.1f uA/%.1f uA) : %.3f ", avg_current_bcm_cut, heep_Ib_simc,  heep_kin1_rates * (avg_current_bcm_cut/heep_Ib_simc) ) << endl;}
	if(abs(6.8-shms_angle)<0.1) { out_file << Form("simc_heep_kin2_rates (shms=6.8 deg) [Hz] x (%.1f uA/%.1f uA) : %.3f ", avg_current_bcm_cut, heep_Ib_simc,  heep_kin2_rates * (avg_current_bcm_cut/heep_Ib_simc) ) << endl;}
	out_file << "" << endl;
	out_file << Form("data_integrated_luminosity [fb^-1]: %.3f", GetLuminosity("data_lumi")) << endl;
	out_file << Form("data_lumiNorm_counts [fb]: %.3f", W_total/GetLuminosity("data_lumi") ) << endl;
	out_file << "" << endl;
	out_file << "# =:=:=:=:=:=:=:=:=:=:=:" << endl;
	out_file << "# SIMC Statistical Goal  " << endl;
	out_file << "# =:=:=:=:=:=:=:=:=:=:=:" << endl;
	out_file << "" << endl;
	if(abs(8.3-shms_angle)<0.1) { 
	  out_file << Form("simc_counts_goal        : %.1f", heep_kin0_counts  ) << endl;
	  out_file << Form("simc_charge_goal [mC]   : %.3f", heep_kin0_charge ) << endl;
	  out_file << Form("simc_integrated_luminosity [fb^-1]: %.4f", GetLuminosity("heep_kin0")) << endl;
	  out_file << Form("simc_lumiNorm_counts [fb]: %.4f", heep_kin0_counts/GetLuminosity("heep_kin0") ) << endl;
	}
	if(abs(7.5-shms_angle)<0.1) { 
	  out_file << Form("simc_counts_goal        : %.1f", heep_kin1_counts  ) << endl;
	  out_file << Form("simc_charge_goal [mC]   : %.3f", heep_kin1_charge ) << endl;
	  out_file << Form("simc_integrated_luminosity [fb^-1]: %.4f", GetLuminosity("heep_kin1")) << endl;
	  out_file << Form("simc_lumiNorm_counts [fb]: %.4f", heep_kin1_counts/GetLuminosity("heep_kin1") ) << endl;
	}
	if(abs(6.8-shms_angle)<0.1) { 
	  out_file << Form("simc_counts_goal        : %.1f", heep_kin2_counts  ) << endl;
	  out_file << Form("simc_charge_goal [mC]   : %.3f", heep_kin2_charge ) << endl;
	  out_file << Form("simc_integrated_luminosity [fb^-1]: %.4f", GetLuminosity("heep_kin2")) << endl;
	  out_file << Form("simc_lumiNorm_counts [fb]: %.4f", heep_kin2_counts/GetLuminosity("heep_kin2") ) << endl;
	}
      }

    if(analysis_cut=="heep_coin")
      {
	out_file << "# =:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:" << endl;
	out_file << "# CaFe H(e,e')p Coincidence Counts  " << endl;
	out_file << "# =:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:" << endl;
    	out_file << "                                     " << endl;
	out_file << Form("heep_total_counts    : %.3f", W_total) << endl;
	out_file << Form("heep_real_counts     : %.3f", W_real)  << endl;
	out_file << Form("heep_random_counts   : %.3f", W_rand)  << endl;
	out_file << "                                     " << endl;
	out_file << Form("heep_real_rate [Hz]  : %.3f", W_real_rate)  << endl;
	out_file << Form("simc_heep_kin0_rates (shms=8.3 deg) [Hz] x (%.1f uA/%.1f uA) : %.3f ", avg_current_bcm_cut, heep_Ib_simc,  heep_kin0_rates * (avg_current_bcm_cut/heep_Ib_simc) ) << endl;	
	out_file << "" << endl;
	out_file << Form("data_integrated_luminosity [fb^-1]: %.3f", GetLuminosity("data_lumi")) << endl;
	out_file << Form("data_lumiNorm_counts [fb]: %.3f", W_real/GetLuminosity("data_lumi") ) << endl;
	out_file << "" << endl;
	out_file << "# =:=:=:=:=:=:=:=:=:=:=:" << endl;
	out_file << "# SIMC Statistical Goal  " << endl;
	out_file << "# =:=:=:=:=:=:=:=:=:=:=:" << endl;
	out_file << Form("simc_counts_goal        : %.1f", heep_kin0_counts  ) << endl;
	out_file << Form("simc_charge_goal [mC]   : %.3f", heep_kin0_charge ) << endl;
	out_file << Form("simc_integrated_luminosity [fb^-1]: %.4f", GetLuminosity("heep_kin0")) << endl;
	out_file << Form("simc_lumiNorm_counts [fb]: %.4f", heep_kin0_counts/GetLuminosity("heep_kin0") ) << endl;
	out_file << "                                     " << endl;
      }
    
    if(analysis_cut=="MF")
      {
	out_file << "# =:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:" << endl;
	out_file << "# CaFe A(e,e')p Mean-Field (MF) Counts  " << endl;
	out_file << "# =:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:" << endl;
	out_file << "                                     " << endl;
	out_file << Form("MF_total_counts    : %.3f", Pm_total) << endl;
	out_file << Form("MF_real_counts     : %.3f", Pm_real )  << endl;
	out_file << Form("MF_random_counts   : %.3f", Pm_rand )  << endl;
	out_file << "                                     " << endl;
	out_file << Form("MF_real_rate [Hz]  : %.3f", Pm_real_rate)  << endl;
	out_file << Form("MF_simc_rate [Hz] x (%.1f uA/%.1f uA) : %.3f ", avg_current_bcm_cut, cafe_Ib_simc,  simc_cafe_rates * (avg_current_bcm_cut/cafe_Ib_simc) ) << endl;	
	out_file << "                                     " << endl;	
	out_file << Form("data_integrated_luminosity [fb^-1]: %.3f", GetLuminosity("data_lumi")) << endl;
	out_file << Form("data_lumiNorm_counts [fb]: %.3f", Pm_real/GetLuminosity("data_lumi") ) << endl;
	out_file << "                                     " << endl;
	out_file << "# =:=:=:=:=:=:=:=:=:=:=:" << endl;
	out_file << "# SIMC Statistical Goal  " << endl;
	out_file << "# =:=:=:=:=:=:=:=:=:=:=:" << endl;
	out_file << Form("simc_counts_goal : %.3f", total_simc_counts ) << endl;
	out_file << Form("simc_charge_goal [mC] : %.3f", total_simc_charge ) << endl;
	out_file << Form("simc_integrated_luminosity [fb^-1]: %.4f", GetLuminosity("simc_lumi")) << endl;
	out_file << Form("simc_lumiNorm_counts [fb]: %.4f", total_simc_counts/GetLuminosity("simc_lumi") ) << endl;
	out_file << "                                     " << endl;


      }
    if(analysis_cut=="SRC")
      {
	out_file << "# =:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:" << endl;
	out_file << "# CaFe A(e,e')p Short-Range Correlated (SRC) Counts  " << endl;
	out_file << "# =:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:" << endl;
	out_file << "                                     " << endl;
	out_file << Form("SRC_total_counts    : %.3f", Pm_total) << endl;
	out_file << Form("SRC_real_counts     : %.3f", Pm_real)  << endl;
	out_file << Form("SRC_random_counts   : %.3f", Pm_rand)  << endl;
	out_file << "                                     " << endl;
	out_file << Form("SRC_real_rate [Hz]  : %.3f", Pm_real_rate)  << endl;
	out_file << Form("SRC_simc_rate [Hz] x (%.1f uA/%.1f uA) : %.3f ", avg_current_bcm_cut, cafe_Ib_simc,  simc_cafe_rates * (avg_current_bcm_cut/cafe_Ib_simc) ) << endl;	
	out_file << "" << endl;
	out_file << Form("data_integrated_luminosity [fb^-1]: %.3f", GetLuminosity("data_lumi")) << endl;
	out_file << Form("data_lumiNorm_counts [fb]: %.3f", Pm_real/GetLuminosity("data_lumi") ) << endl;
	out_file << "" << endl;
	out_file << "# =:=:=:=:=:=:=:=:=:=:=:" << endl;
	out_file << "# SIMC Statistical Goal  " << endl;
	out_file << "# =:=:=:=:=:=:=:=:=:=:=:" << endl;
	out_file << Form("simc_counts_goal : %.3f", total_simc_counts ) << endl;
	out_file << Form("simc_charge_goal [mC] : %.3f", total_simc_charge ) << endl;
	out_file << Form("simc_integrated_luminosity [fb^-1]: %.4f", GetLuminosity("simc_lumi")) << endl;
	out_file << Form("simc_lumiNorm_counts [fb]: %.4f", total_simc_counts/GetLuminosity("simc_lumi") ) << endl;
	out_file << "                                     " << endl;
      }

    if(analysis_cut!="bcm_calib"){
    out_file << "                                     " << endl;
    out_file << "# =:=:=:=:=:=:=:=:=:=:=:=:" << endl;
    out_file << "# DAQ Trigger Information  " << endl;
    out_file << "# =:=:=:=:=:=:=:=:=:=:=:=:" << endl;
    out_file << "                                     " << endl;
    out_file << "# NOTE: scaler triggers are not pre-scaled | accepted triggers are pre-scaled  " << endl;
    out_file << "#       cpu_live_time   = T#_accepted / ( T#_scaler / Ps#_factor)       " << endl;
    out_file << "#       total_live_time = edtm_accepted / ( edtm_scaler / Ps#_factor)       " << endl;
    out_file << "                                     " << endl;
    out_file << Form("edtm_scaler   :  %.3f  ",         total_edtm_scaler_bcm_cut ) << endl;
    out_file << Form("edtm_accepted :  %.3f  ",         total_edtm_accp_bcm_cut) << endl;
    out_file << "                                     " << endl;
    out_file << "# pre-scale factors (-1: trigger OFF)                 " << endl;
    out_file << Form("Ps1_factor: %.1f", Ps1_factor) << endl;
    out_file << Form("Ps2_factor: %.1f", Ps2_factor) << endl;
    out_file << Form("Ps3_factor: %.1f", Ps3_factor) << endl;
    out_file << Form("Ps4_factor: %.1f", Ps4_factor) << endl;
    out_file << Form("Ps5_factor: %.1f", Ps5_factor) << endl;
    out_file << Form("Ps6_factor: %.1f", Ps6_factor) << endl;
    out_file << "                                     " << endl;
    out_file << "# pre-trigger scalers                 " << endl;
    out_file << Form("S1X_scaler:  %.3f [ %.3f kHz ] ", total_s1x_scaler_bcm_cut,    S1XscalerRate_bcm_cut) << endl; 
    out_file << Form("T1_scaler:  %.3f [ %.3f kHz ] ",  total_trig1_scaler_bcm_cut,  TRIG1scalerRate_bcm_cut) << endl;
    out_file << Form("T2_scaler:  %.3f [ %.3f kHz ] ",  total_trig2_scaler_bcm_cut,  TRIG2scalerRate_bcm_cut) << endl;
    out_file << Form("T3_scaler:  %.3f [ %.3f kHz ] ",  total_trig3_scaler_bcm_cut,  TRIG3scalerRate_bcm_cut) << endl;
    out_file << Form("T4_scaler:  %.3f [ %.3f kHz ] ",  total_trig4_scaler_bcm_cut,  TRIG4scalerRate_bcm_cut) << endl;
    out_file << Form("T5_scaler:  %.3f [ %.3f kHz ] ",  total_trig5_scaler_bcm_cut,  TRIG5scalerRate_bcm_cut) << endl;
    out_file << Form("T6_scaler:  %.3f [ %.3f kHz ] ",  total_trig6_scaler_bcm_cut,  TRIG6scalerRate_bcm_cut) << endl;
    
    out_file << "                                     " << endl;
    out_file << "# accepted triggers (pre-scaled)     " << endl;
    out_file << Form("T1_accepted: %.3f [ %.3f kHz ]  ", total_trig1_accp_bcm_cut,    TRIG1accpRate_bcm_cut) << endl;
    out_file << Form("T2_accepted: %.3f [ %.3f kHz ]  ", total_trig2_accp_bcm_cut,    TRIG2accpRate_bcm_cut) << endl;
    out_file << Form("T3_accepted: %.3f [ %.3f kHz ]  ", total_trig3_accp_bcm_cut,    TRIG3accpRate_bcm_cut) << endl;
    out_file << Form("T4_accepted: %.3f [ %.3f kHz ]  ", total_trig4_accp_bcm_cut,    TRIG4accpRate_bcm_cut) << endl;
    out_file << Form("T5_accepted: %.3f [ %.3f kHz ]  ", total_trig5_accp_bcm_cut,    TRIG5accpRate_bcm_cut) << endl;
    out_file << Form("T6_accepted: %.3f [ %.3f kHz ]  ", total_trig6_accp_bcm_cut,    TRIG6accpRate_bcm_cut) << endl;
    out_file << "                                     " << endl;
    out_file << "# daq computer (cpu) and total live time  " << endl;
    if(Ps1_factor > -1) {
      out_file << Form("T1_cpuLT:    %.3f +- %.3f ",  cpuLT_trig1,                 cpuLT_trig1_err_Bi) << endl;
      out_file << Form("T1_tLT:      %.3f +- %.3f ",  tLT_trig1,                   tLT_trig1_err_Bi) << endl;	
      out_file << "                                     " << endl;
    }
    if(Ps2_factor > -1) {
      out_file << Form("T2_cpuLT:    %.3f +- %.3f ",  cpuLT_trig2,                 cpuLT_trig2_err_Bi) << endl;
      out_file << Form("T2_tLT:      %.3f +- %.3f ",  tLT_trig2,                   tLT_trig2_err_Bi) << endl;	
      out_file << "                                     " << endl;
    }
    if(Ps3_factor > -1) {
      out_file << Form("T3_cpuLT:    %.3f +- %.3f ",  cpuLT_trig3,                 cpuLT_trig3_err_Bi) << endl;
      out_file << Form("T3_tLT:      %.3f +- %.3f ",  tLT_trig3,                   tLT_trig3_err_Bi) << endl;	
      out_file << "                                     " << endl;
    }
    if(Ps4_factor > -1) {
      out_file << Form("T4_cpuLT:    %.3f +- %.3f ",  cpuLT_trig4,                 cpuLT_trig4_err_Bi) << endl;
      out_file << Form("T4_tLT:      %.3f +- %.3f ",  tLT_trig4,                   tLT_trig4_err_Bi) << endl;	
      out_file << "                                     " << endl;
    }
    if(Ps5_factor > -1) {
      out_file << Form("T5_cpuLT:    %.3f +- %.3f ",  cpuLT_trig5,                 cpuLT_trig5_err_Bi) << endl;
      out_file << Form("T5_tLT:      %.3f +- %.3f ",  tLT_trig5,                   tLT_trig5_err_Bi) << endl;	
      out_file << "                                     " << endl;
    }
    if(Ps6_factor > -1) {
      out_file << Form("T6_cpuLT:    %.3f +- %.3f ",  cpuLT_trig6,                 cpuLT_trig6_err_Bi) << endl;
      out_file << Form("T6_tLT:      %.3f +- %.3f ",  tLT_trig6,                   tLT_trig6_err_Bi) << endl;	
      out_file << "                                     " << endl;
    }
    out_file << "# =:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:" << endl;
    out_file << "# Drift Chambers Tracking Efficiency  " << endl;
    out_file << "# =:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:" << endl;
    out_file << "                                     " << endl;
    out_file << Form("hms_had_track_eff:  %.3f +- %.3f",  hTrkEff,  hTrkEff_err) << endl;
    out_file << Form("shms_elec_track_eff: %.3f +- %.3f",  pTrkEff, pTrkEff_err) << endl;
    out_file << "                                     " << endl;
    out_file << "# =:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:" << endl;
    out_file << "# Data Analysis Cuts                  " << endl;
    out_file << "# =:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:" << endl;
    out_file << "                                   " << endl;
    out_file << "#--- Tracking Efficiency Definition --- " << endl;
    out_file << "# tracking efficiency = (did && should) / should" << endl;
    out_file << "" << endl;
    if((analysis_cut=="heep_coin") || (analysis_cut=="MF") || (analysis_cut=="SRC") )
      {
	if(hdc_ntrk_cut_flag)    {out_file << Form("# (did) HMS min. number of tracks (H.dc.ntrack): >= %.1f", c_hdc_ntrk_min) << endl;}
	if(hScinGood_cut_flag)   {out_file <<      "# (should) HMS good (fiducial) scintillator hit (H.hod.goodscinhit): true"  << endl;}
	if(hcer_cut_flag)        {out_file << Form("# (should) HMS gas Cherenkov number of photoelectrons (H.cer.npeSum): (%.1f, %.1f)", c_hnpeSum_min, c_hnpeSum_max) << endl;}
	if(hetotnorm_cut_flag)   {out_file << Form("# (should) HMS calorimeter energy / central_momentum  (H.cal.etotnorm): (%.1f, %.1f)", c_hetotnorm_min, c_hetotnorm_max) << endl;}
	if(hBeta_notrk_cut_flag) {out_file << Form("# (should) HMS hodoscope beta no_track (H.hod.betanotrack): (%.1f. %.1f)", c_hBetaNtrk_min, c_hBetaNtrk_max) << endl;}
      }
    if((analysis_cut=="heep_singles") || (analysis_cut=="heep_coin") || (analysis_cut=="MF") || (analysis_cut=="SRC") )
      {
	out_file << "                                   " << endl;
	if(pdc_ntrk_cut_flag)    {out_file << Form("# (did) SHMS min. number of tracks (P.dc.ntrack): >= %.1f", c_pdc_ntrk_min) << endl;}
	if(pScinGood_cut_flag)   {out_file <<      "# (should) SHMS good (fiducial) scintillator hit (P.hod.goodscinhit): true"  << endl;}
	if(pngcer_cut_flag)      {out_file << Form("# (should) SHMS noble gas Chrenkov number of photoelectrons (P.ngcer.npeSum): (%.1f, %.1f)", c_pngcer_npeSum_min, c_pngcer_npeSum_max) << endl;}
	if(phgcer_cut_flag)      {out_file << Form("# (should) SHMS heavy gas Chrenkov number of photoelectrons (P.hgcer.npeSum): (%.1f, %.1f)", c_phgcer_npeSum_min, c_phgcer_npeSum_max) << endl;}
	if(petotnorm_cut_flag)   {out_file << Form("# (should) SHMS calorimeter energy / central_momentum  (p.cal.etotnorm): (%.1f, %.1f)", c_petotnorm_min, c_petotnorm_max) << endl;}
	if(pBeta_notrk_cut_flag) {out_file << Form("# (should) SHMS hodoscope beta no_track (P.hod.betanotrack): (%.1f. %.1f)", c_pBetaNtrk_min, c_pBetaNtrk_max) << endl;}
      }
    out_file << "                                   " << endl;
      if(ePctime_cut_flag && ((analysis_cut=="heep_coin") || (analysis_cut=="MF") || (analysis_cut=="SRC") ))     {
      out_file << "#---Coincidence Time Cut--- " << endl;
      out_file << Form("# electron (SHMS)-proton(HMS) (prompt) coincidence time (CTime.epCoinTime_ROC2):   (%.3f, %.3f) [ns]", ePctime_cut_min, ePctime_cut_max) << endl;
      out_file << Form("# electron (SHMS)-proton(HMS) (left)   accidentals sample: (%.3f, %.3f) [ns]", ePctime_cut_max_L, ePctime_cut_min_L) << endl;
      out_file << Form("# electron (SHMS)-proton(HMS) (right)  accidentals sample: (%.3f, %.3f) [ns]", ePctime_cut_min_R, ePctime_cut_max_R) << endl;      
    }
    out_file << "                                   " << endl;
    out_file << "#---Acceptance Cuts--- " << endl;
    if((analysis_cut=="heep_coin")  ||  (analysis_cut=="MF") || (analysis_cut=="SRC"))
       {    
	 if(hdelta_cut_flag)        {out_file << Form("# HMS Momentum Acceptance (H.gtr.dp): (%.3f, %.3f) [%%]",  c_hdelta_min, c_hdelta_max ) << endl;}
	 if(hxptar_cut_flag)        {out_file << Form("# HMS Out-of-Plane (xptar) Angular Acceptance (H.gtr.th): (%.3f, %.3f) [radians]",  c_hxptar_min, c_hxptar_max ) << endl;}
	 if(hyptar_cut_flag)        {out_file << Form("# HMS In-Plane (yptar) Angular Acceptance (H.gtr.ph): (%.3f, %.3f) [radians]",  c_hyptar_min, c_hyptar_max ) << endl;}
	 if(hmsCollCut_flag)        {out_file << "# HMS Collimator Cut: ON " << endl;}
	 
	 if(edelta_cut_flag)        {out_file << Form("# SHMS Momentum Acceptance (P.gtr.dp): (%.3f, %.3f) [%%]", c_edelta_min, c_edelta_max ) << endl;}
	 if(exptar_cut_flag)        {out_file << Form("# SHMS Out-of-Plane (xptar) Angular Acceptance (P.gtr.th): (%.3f, %.3f) [radians]",  c_exptar_min, c_exptar_max ) << endl;}
	 if(eyptar_cut_flag)        {out_file << Form("# SHMS In-Plane (yptar) Angular Acceptance (P.gtr.ph): (%.3f, %.3f) [radians]",  c_eyptar_min, c_eyptar_max ) << endl;}
	 if(shmsCollCut_flag)       {out_file << "# SHMS Collimator Cut: ON " << endl;}
	
	 if(ztarDiff_cut_flag)      {out_file << Form("# Z-Reaction Vertex Difference (H.react.z-P.react.z): (%.3f, %.3f) [cm]", c_ztarDiff_min, c_ztarDiff_max ) << endl;}
	
       }
    if((analysis_cut=="heep_singles") || (analysis_cut=="optics")  ||  (analysis_cut=="lumi") )
      {
	if(edelta_cut_flag)        {out_file << Form("# SHMS Momentum Acceptance (P.gtr.dp): (%.3f, %.3f) [%%]", c_edelta_min, c_edelta_max ) << endl;}
	if(exptar_cut_flag)        {out_file << Form("# SHMS Out-of-Plane (xptar) Angular Acceptance (P.gtr.th): (%.3f, %.3f) [radians]",  c_exptar_min, c_exptar_max ) << endl;}
	if(eyptar_cut_flag)        {out_file << Form("# SHMS In-Plane (yptar) Angular Acceptance (P.gtr.ph): (%.3f, %.3f) [radians]",  c_eyptar_min, c_eyptar_max ) << endl;}
	if(shmsCollCut_flag)       {out_file << "# SHMS Collimator Cut: ON " << endl;}


      }
    out_file << "#                       " << endl;
    out_file << "#---Particle Identification (PID) Cuts--- " << endl;
    if((analysis_cut=="heep_coin")  ||  (analysis_cut=="MF") || (analysis_cut=="SRC"))
      {    
	if(hetot_trkNorm_pidCut_flag) {out_file << Form("# HMS calorimeter total energy / track momentum (H.cal.etottracknorm): (%.1f, %.1f)", cpid_hetot_trkNorm_min, cpid_hetot_trkNorm_max) << endl;}
	if(hcer_pidCut_flag)          {out_file << Form("# HMS gas Cherenkov number of photoelectrons (H.cer.npeSum): (%.1f, %.1f)", cpid_hcer_npeSum_min, cpid_hcer_npeSum_max) << endl;}
      }
      if((analysis_cut=="heep_singles") || (analysis_cut=="heep_coin") ||  (analysis_cut=="MF") || (analysis_cut=="SRC"))
      {
	if(petot_trkNorm_pidCut_flag)  {out_file << Form("# SHMS calorimeter total energy / track momentum (P.cal.etottracknorm): (%.1f, %.1f)",cpid_petot_trkNorm_min, cpid_petot_trkNorm_max) << endl;}
	if(pngcer_pidCut_flag)        {out_file << Form("# SHMS noble gas Chrenkov number of photoelectrons (P.ngcer.npeSum): (%.1f, %.1f)",cpid_pngcer_npeSum_min,cpid_pngcer_npeSum_max) << endl;}
	if(phgcer_pidCut_flag)        {out_file << Form("# SHMS heavy gas Chrenkov number of photoelectrons (P.hgcer.npeSum): (%.1f, %.1f)",cpid_phgcer_npeSum_min,cpid_phgcer_npeSum_max) << endl;}
      }
    out_file << "#                       " << endl;
    out_file << "#---Kinematics Cuts--- " << endl;
    if((analysis_cut=="heep_singles") || (analysis_cut=="heep_coin"))
      {
	if(Q2_heep_cut_flag)  {out_file << Form("# H(e,e'p) 4-momentum transferred squared, Q2 (P.kin.primary.Q2): (%.3f, %.3f) [GeV2]", c_heep_Q2_min, c_heep_Q2_max) << endl;}
	if(xbj_heep_cut_flag) {out_file << Form("# H(e,e'p) x-Bjorken, Xbj (P.kin.primary.x_bj): (%.3f, %.3f)", c_heep_xbj_min, c_heep_xbj_max) << endl;}
	if(W_heep_cut_flag)  {out_file << Form("# H(e,e'p) Invariant Mass, W (P.kin.primary.W): (%.3f, %.3f) [GeV]", c_heep_W_min, c_heep_W_max) << endl;}
      }
     if(analysis_cut=="heep_coin")
       {
	 if(MM_heep_cut_flag)  {out_file << Form("# H(e,e'p) Missing Mass, MM (H.kin.secondary.Mrecoil): (%.3f, %.3f) [GeV]", c_heep_MM_min, c_heep_MM_max) << endl;}
	 if(Em_heep_cut_flag)  {out_file << Form("# H(e,e'p) Missing Energy, Em=nu-Ep (H.kin.secondary.Em): (%.3f, %.3f) [GeV]", c_heep_Em_min, c_heep_Em_max) << endl;}
       }
     if(analysis_cut=="MF")
       {
	 if(Q2_MF_cut_flag) {out_file << Form("# A(e,e'p) 4-momentum transferred squared, Q2 (P.kin.primary.Q2): (%.3f, %.3f) [GeV2]", c_MF_Q2_min, c_MF_Q2_max) << endl;}
	 if(Pm_MF_cut_flag) {out_file << Form("# A(e,e'p) Missing Momentum, Pm (H.kin.secondary.Pm): (%.3f, %.3f) [GeV]", c_MF_Pm_min, c_MF_Pm_max) << endl;}
	 if(Em_d2MF_cut_flag && tgt_type=="LD2") {out_file << Form("# A(e,e'p) Missing Energy, Em (H.kin.secondary.Em_nuc): (%.3f, %.3f) [GeV]", c_d2MF_Em_min, c_d2MF_Em_max) << endl;}
	 if(Em_MF_cut_flag && tgt_type!="LD2") {out_file << Form("# A(e,e'p) Missing Energy, Em (H.kin.secondary.Em_nuc): (%.3f, %.3f) [GeV]", c_MF_Em_min, c_MF_Em_max) << endl;}

       }
     if(analysis_cut=="SRC")
       {
	 if(Q2_SRC_cut_flag) {out_file << Form("# A(e,e'p) 4-momentum transferred squared, Q2 (P.kin.primary.Q2): (%.3f, %.3f) [GeV2]", c_SRC_Q2_min, c_SRC_Q2_max) << endl;}
	 if(Pm_SRC_cut_flag) {out_file << Form("# A(e,e'p) Missing Momentum, Pm (H.kin.secondary.Pm): (%.3f, %.3f) [GeV]", c_SRC_Pm_min, c_SRC_Pm_max) << endl;}
	 if(Xbj_SRC_cut_flag) {out_file << Form("# A(e,e'p) x-Bjorken, Xbj (P.kin.primary.x_bj): (%.3f, %.3f)", c_SRC_Xbj_min, c_SRC_Xbj_max) << endl;}
	 if(thrq_SRC_cut_flag) {out_file << Form("# A(e,e'p) theta_rq  (H.kin.secondary.th_bq): (%.3f, %.3f) [deg]", c_SRC_thrq_min, c_SRC_thrq_max) << endl;}	 
	 if(Em_d2SRC_cut_flag && tgt_type=="LD2") {out_file << Form("# A(e,e'p) Missing Energy, Em (H.kin.secondary.Em_nuc): (%.3f, %.3f) [GeV]", c_d2MF_Em_min, c_d2MF_Em_max) << endl;}
	 if(Em_SRC_cut_flag && tgt_type!="LD2") {out_file << "# A(e,e'p) Dynamic Missing Energy (A>2 nuclei), Em_src = nu - Tp - (sqrt(MN*MN + Pm*Pm) - MN) | see definition in baseAnalyzer.cpp" << endl;}

       }
					    
    } // end !bcm_calib requirement
    
    
    // CLOSE files
    out_file.close();
    in_file.close();
    
  }
  
}

//_______________________________________________________________________________
void baseAnalyzer::WriteReportSummary()
{
  
  /*Method to write charge, efficiencies, live time and other relevant quantities to a data file
    on a run-by-run basis, and self-updating file, meaning, each run that is replayed will be appended into the file.    
   */
  
  cout << "Calling WriteReportSummary() . . ." << endl;

  
  if(analyze_data==true){

    //---------------------------------------------------------

    
    //Check if file already exists
    in_file.open(output_SummaryFileName.Data());
    
    if(in_file.fail()){
      
      cout << "Report File does NOT exist, will create one . . . " << endl;
      
      out_file.open(output_SummaryFileName);
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
      if(ePctime_cut_flag)         {out_file << Form("# Proton Coincidence Time Cut: (%.3f, %.3f )ns", ePctime_cut_min, ePctime_cut_max) << endl;}
      if(petot_trkNorm_pidCut_flag)   {out_file << Form("# SHMS Calorimeter EtotTrackNorm Cut: (%.3f, %.3f)", cpid_petot_trkNorm_min,  cpid_petot_trkNorm_max) << endl;}
      if(pngcer_pidCut_flag) {out_file << Form("# SHMS Noble Gas Cherenkov NPE Sum Cut: (%.3f, %.3f)", cpid_pngcer_npeSum_min,  cpid_pngcer_npeSum_max) << endl;}
      if(phgcer_pidCut_flag) {out_file << Form("# SHMS Heavy Gas Cherenkov NPE Sum Cut: (%.3f, %.3f)", cpid_phgcer_npeSum_min,  cpid_phgcer_npeSum_max) << endl;}
      if(hetot_trkNorm_pidCut_flag) {out_file << Form("# HMS Calorimeter EtotTrackNorm Cut: (%.3f, %.3f)", cpid_hetot_trkNorm_min,  cpid_hetot_trkNorm_max) << endl;}
      if(hcer_pidCut_flag) {out_file << Form("# HMS Gas Cherenkov NPE Sum Cut: (%.3f, %.3f)", cpid_hcer_npeSum_min,  cpid_hcer_npeSum_max) << endl;}      
      out_file << "#                                     " << endl;
      out_file << "#---Kinematics Cuts--- " << endl;
      if(Q2_heep_cut_flag)            {out_file << Form("# 4-Momentum Transfer (Q^2): (%.3f, %.3f) GeV^2", c_heep_Q2_min, c_heep_Q2_max ) << endl;}
      if(Em_heep_cut_flag)            {out_file << Form("# Missing Energy, Em: (%.3f, %.3f) GeV",   c_heep_Em_min, c_heep_Em_max ) << endl;}
      if(W_heep_cut_flag)             {out_file << Form("# Invariant Mass, W: (%.3f, %.3f) GeV",   c_heep_W_min,  c_heep_W_max  ) << endl;}
      if(MM_heep_cut_flag)          {out_file << Form("# Missing Mass, MM: (%.3f, %.3f) GeV",   c_heep_MM_min,  c_heep_MM_max  ) << endl;}
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
    out_file.open(output_SummaryFileName, ios::out | ios::app);
    out_file << std::setw(7) << run  << std::setw(25) << total_charge_bcm_cut << std::setw(25) << avg_current_bcm_cut << std::setw(25) << hTrkEff << std::setw(25) << hTrkEff_err << std::setw(25) << pTrkEff << std::setw(25) << pTrkEff_err << std::setw(25) << tgtBoil_corr << std::setw(25) << tgtBoil_corr_err << std::setw(25) << hadAbs_corr << std::setw(25) << hadAbs_corr_err << std::setw(25) << cpuLT_trig << std::setw(25) << cpuLT_trig_err_Bi << std::setw(25) << cpuLT_trig_err_Bay << std::setw(25) << tLT_trig << std::setw(25) << tLT_trig_err_Bi << std::setw(25) << tLT_trig_err_Bay << std::setw(25) << S1XscalerRate_bcm_cut << std::setw(25) << trig_rate << std::setw(25) << EDTMscalerRate_bcm_cut << std::setw(25) << Ps_factor << std::setw(25) << total_edtm_accp_bcm_cut << std::setw(25) << (total_edtm_scaler_bcm_cut / Ps_factor) << std::setw(25) << total_trig_accp_bcm_cut << std::setw(25) << (total_trig_scaler_bcm_cut / Ps_factor) << endl;
    out_file.close();
  }
  
  cout << "Ending WriteReportSummary() . . ." << endl;
  
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

//______________________________________________________________________________
Double_t baseAnalyzer::GetLuminosity(TString user_input="")
{
  /* 
     Brief: calculates luminosity  as follows: luminosity = Constant * total_charge / targetfac,
     where targetfac = m_amu / tgt_areal_density
     
     Standard definition:  differential luminosity, dL = (1/xsec) * dN/dt, where xsec is the cross section 
     and dN/dt is the event count rate. If 'xsec' is in [cm^2] and 'dN/dt' is in [Hz], then 
     luminosity is in [cm^-2 sec^-1]. One can also get the integrated luminosity by integrated over
     total beam-on-target run time.

     To get the constant in units of cm^2, given the areal density in g/cm^2, and charge in milliCoulomb :
     L = Constant * (1 mC / e-) *  (1 g / cm^2) / 1 amu
     
     ** 1 e- = 1.60217663 × 10^-19 Coulombs ---> 1 mC / e- = 1/ (1.60217663 × 10^-19 ) * 1e-3 mC/ 1 C = 6241509090043338.0 ---> 1 mC/e- ~ 6.24e15
     
     ** (1 g / cm^2) / 1 amu  =  NA atoms / mol = 6.02214076 x 10^23

     C = (6.241509090043338 x 10^15)  * (6.02214076 x 10^23) = 3.7587246295060495e+39 cm^-2  | 1 ubarn = 1e-30 cm^2
     
     C[ub^-1] =  3.7587246295060495e+39 cm^-2 * 1e-30 cm^2 / 1 ub = 3758724629.5060496 ----> C ~ 3.75872 x 10^9 ub^-1  (inverse microbarns)
     
      might be more useful to convert to femtobarns, just in case: here it is
     C[cm^-2] * 1 cm^2 / 1e39 [fb]

     or it might be better to convert from cm^-2 to  GeV^2
     C[cm^-2] * (1 cm^2 / 1e36 [pb])  * (1 [pb] / 2.56819×10−9 GeV^−2) = 1. [cm^2] / 3.8937929e-28 [GeV^-2] ;

     To directly compare different targets, one can normalize the total experimental counts (N)  by total luminosity

   */
  
  // initialize return values
  luminosity =  luminosity_simc = heep_kin0_lumi_simc = heep_kin1_lumi_simc = heep_kin2_lumi_simc = 0.0;

  
  // [g/cm^2]     =    [g/cm^3]    *   [cm]
  tgt_areal_density =  tgt_density *  tgt_thickness;  

  //  [cm^-2]      [g/cm^2]      /  [g/mol]  --> becomaes unitless, since constant was extracted from this
  targetfac = tgt_areal_density / tgt_mass ;  

  Double_t Constant = NA * (1./elementary_charge) * 1e-3 * 1e-39; // units: [fb^-1] inverse-femtobarns
  
  // calculate DATA integrated luminosity
  //              [fb^-1]        [mC]             [unitless]
  luminosity =  Constant * total_charge_bcm_cut * targetfac;  

  if( (analysis_cut=="MF") || (analysis_cut=="SRC") ){

    // calculate simc integrated luminosity (based on simulation predictions)
    total_simc_time = stod(split(FindString(Form("%s_%s_time", tgt_type.Data(), analysis_cut.Data()),    input_SIMCinfo_FileName.Data())[0], '=')[1]); // [hr]
    cafe_Ib_simc = stod(split(FindString("cafe_Ib_simc",    input_SIMCinfo_FileName.Data())[0], '=')[1]); //[uA]
    
    // [mC]                [uC / sec]        [hr]      [sec]/[hr]  0.001 mC / 1 uC
    total_simc_charge =  cafe_Ib_simc * total_simc_time * 3600. * 1e-3;
    
    luminosity_simc = Constant * total_simc_charge * targetfac;

  }
  else if ((analysis_cut=="heep_singles") || (analysis_cut=="heep_coin")){

    heep_kin0_time = stod(split(FindString("heep_kin0_time",    input_SIMCinfo_FileName.Data())[0], '=')[1]); // [hr]
    heep_kin1_time = stod(split(FindString("heep_kin1_time",    input_SIMCinfo_FileName.Data())[0], '=')[1]); // [hr]
    heep_kin2_time = stod(split(FindString("heep_kin2_time",    input_SIMCinfo_FileName.Data())[0], '=')[1]); // [hr]

    heep_Ib_simc = stod(split(FindString("heep_Ib_simc",    input_SIMCinfo_FileName.Data())[0], '=')[1]); //[uA]

    // [mC]                [uC / sec]        [hr]      [sec]/[hr]  0.001 mC / 1 uC
    heep_kin0_charge =  heep_Ib_simc * heep_kin0_time * 3600. * 1e-3;
    heep_kin1_charge =  heep_Ib_simc * heep_kin1_time * 3600. * 1e-3;
    heep_kin2_charge =  heep_Ib_simc * heep_kin2_time * 3600. * 1e-3;

    heep_kin0_lumi_simc = Constant * heep_kin0_charge * targetfac;
    heep_kin1_lumi_simc = Constant * heep_kin1_charge * targetfac;
    heep_kin2_lumi_simc = Constant * heep_kin2_charge * targetfac;
    
  }
  

  if(user_input=="data_lumi"){
    return luminosity;  // [fb^-1] 
  }
  else if(user_input=="simc_lumi"){
    return luminosity_simc;  // [fb^-1]     
  }
  else if(user_input=="heep_kin0"){
    return heep_kin0_lumi_simc; // [fb^-1]
  }
  else if(user_input=="heep_kin1"){
    return heep_kin1_lumi_simc; // [fb^-1]
  }
  else if(user_input=="heep_kin2"){
    return heep_kin2_lumi_simc; // [fb^-1]
  }
  else{
    return 0.0;
  }
  
}

//______________________________________________________________________________
void baseAnalyzer::MakePlots()
{
  cout << "Calling MakePlots() . . . " << endl;

  string cmd0 = Form("emacs -nw %s", output_ReportFileName.Data());
  cout << cmd0.c_str() << endl;
  gSystem->Exec(cmd0.c_str());
  
  string cmd=Form("root -l -q -b \"UTILS_CAFE/online_scripts/make_online_plots.cpp(%d, \\\"%s\\\", \\\"%s\\\", \\\"%s\\\", \\\"%s\\\")\" ", run, tgt_type.Data(), analysis_type.Data(), analysis_cut.Data(), data_OutputFileName.Data());
  cout << cmd.c_str() << endl;

  if(analysis_cut!="optics"){
    gSystem->Exec(cmd.c_str());
  }
  
}

//--------------------------MAIN ANALYSIS FUNCTIONS-----------------------------
void baseAnalyzer::run_data_analysis()
{
  /*
    Brief:  This method call all the necessary methods to carry out the 
    full data analysis. 

    This is supposed to be a generic baseAnalyzer class which analyzes data in
    a generic way. The analyzer assumes that all the calibrations from the data
    have been done. See each of the methods for details of what the analyzer does.

    Additional methods may be added in accordance with the necessity of the experiment.
    For example, methods to apply radiative corrections and get cross section still need 
    to be added.

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
  //WriteReportSummary();
  CombineHistos();
  MakePlots();
  
  //------------------

  
}

//--------------------------MAIN ANALYSIS FUNCTIONS-----------------------------
void baseAnalyzer::run_cafe_scalers()
{
 
  //------------------
  ReadInputFile();
  ReadReport();
  
  ReadScalerTree();   
  ScalerEventLoop();       
  CalcEff();
  WriteReport();

  //------------------

  
}
