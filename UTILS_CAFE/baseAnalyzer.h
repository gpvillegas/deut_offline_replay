/*
Author: Carlos Yero
email: cyero002@fiu.edu, cyero@jlab.org
Date Created: August 22, 2020
*/

#ifndef BASE_ANALYZER_H
#define BASE_ANALYZER_H

#include "./UTILS/parse_utils.h" //useful C++ string parsing utilities
#include "./UTILS/hist_utils.h" //useful C++ histogram bin extraction utility

class baseAnalyzer
{
  
public:
  
  //Constructor / Destructor
  baseAnalyzer( int irun=-1, int ievt=-1, string mode="", string earm="", Bool_t ana_data=0, string ana_cuts="", string ana_type="", Bool_t hel_flag=0, string bcm_name="", double thrs=-1, string trig="", Bool_t combine_flag=0); //initialize member variables
  ~baseAnalyzer();
  
  //MAIN ANALYSIS FUNCTIONS
  void run_data_analysis();
  void run_cafe_scalers(); // mainly for generating cafe output file (for bcm calib runs)
  
  //Function prototypes
  void ReadInputFile();
  void ReadReport();
  void SetHistBins();
  void CreateHist();
  void ReadScalerTree();  
  void ScalerEventLoop(); //bcm current cut threshold in uA units
  void ReadTree();
  void EventLoop();
  void CalcEff();
  void ApplyWeight();
  void WriteHist();
  void WriteReport();
  void WriteReportSummary();
  void CombineHistos();
  
  //void CalcRadCorr(); 
  //void ApplyRadCorr();
  //void ChargeNorm(); 
  void RandSub(); //Apply subtraction of random coincidence background
  //void GetAsymmetry();
  
  // Helper Functions
  Double_t GetCoinTimePeak();
  void MakePlots();
  
protected:

  //Set Constants
  const Double_t pi = TMath::Pi(); 
  const Double_t dtr = pi / 180.;
  const Double_t amu2GeV = 0.93149432;

  // target mass (amu) 
  Double_t MH_amu     = 1.00794       ;
  Double_t MD_amu     = 2.01410177812 ;
  Double_t MBe9_amu   = 9.012182      ;
  Double_t MB10_amu   = 10.0129370    ;
  Double_t MB11_amu   = 11.009306     ;
  Double_t MC12_amu   = 12.0107       ;
  Double_t MAl27_amu  = 26.98153      ;
  Double_t MCa40_amu  = 39.962590863  ;
  Double_t MCa48_amu  = 47.95252276   ;
  Double_t MFe54_amu  = 53.9396147    ;
  Double_t MTi48_amu  = 47.9479463    ;
  // target mass (GeV)
  Double_t MH     = MH_amu    * amu2GeV;
  Double_t MD     = MD_amu    * amu2GeV;
  Double_t MBe9   = MBe9_amu  * amu2GeV;
  Double_t MB10   = MB10_amu  * amu2GeV;
  Double_t MB11   = MB11_amu  * amu2GeV;
  Double_t MC12   = MC12_amu  * amu2GeV;
  Double_t MAl27  = MAl27_amu * amu2GeV;
  Double_t MCa40  = MCa40_amu * amu2GeV;
  Double_t MCa48  = MCa48_amu * amu2GeV;
  Double_t MFe54  = MFe54_amu * amu2GeV;
  Double_t MTi48  = MTi48_amu * amu2GeV;
  
  //detected particle masses (GeV/c^2)
  const Double_t me = 0.000510998950;  //electron Mass
  const Double_t MP = 0.938272;  //Proton Mass
  const Double_t MN = 0.939565;  //Neutron Mass
  const Double_t MK = 0.493677;  //Kaon (K+, K-) Mass 
  const Double_t MPi = 0.139570; //Pion (Pi+, Pi-) Mass
  const Double_t MLam = 1.115683; //Neutral Lambda (uds) Mass
  const Double_t MSig = 1.192642; //Neutral Sigma (uds) Mass
  
  
  //Initialization parameters (variables actually used in baseAnalyzer.cpp)
  int run;          // run number
  int evtNum;       // number of events replayed (by the input ROOTfile)
  TString daq_mode;   //"coin" or "singles"
  TString e_arm_name;   // electron arm: "HMS" or "SHMS"
  TString h_arm_name;   //hadron arm
  Bool_t analyze_data;    // analyze data ? if true (analyze data), if false, (analyze_simc)
  TString analysis_cut;    // analysis cuts: either "bcm_calib", "lumi", "optics", "heep_singles", "heep_coin", "MF" or "SRC" 
  TString analysis_type;   // online analysis prefix: either "prod" (for production) or "sample" for sample replay of 100k events (or any other sample evts)
  Bool_t helicity_flag;     //helicity flag
  TString bcm_type;       // BCM type : "BCM1, BCM2, BCM4A, BCM4B, BCM4C"
  Double_t bcm_thrs;      // BCM current threshold cut (analyze data and scalers ONLY above a certain bcm_thrs, e.g. > 5 uA)
  TString trig_type;      // trigger type to actually use when calculating live tim
  Bool_t combine_runs_flag;     //flag to combine multiple runs (usually sequential runs @ same kinematics in an experiment)

  // Read in general info from REPORT file 

  // target type (will be read from report file, rather than user input -- SAFER THIS WAY! :) )
  TString tgt_type;
  Double_t tgt_mass;
  
  Double_t beam_energy;
  Double_t hms_part_mass;
  Double_t hms_p;
  Double_t hms_angle;
  Double_t shms_part_mass;
  Double_t shms_p;
  Double_t shms_angle;
  Double_t run_len;
  
  //Spectrometer prefixes to be used in SetBranchAddress()
  TString eArm;
  TString hArm;
  TString e_arm;
  TString h_arm;
  TString nroc;
  TString daq; //used for storing either "shms", "hms" or "coin"
  TString scl_tree_name;

  //Option to combine multiple runs
  Bool_t combine_histos;
  
  //Declare TFile Pointers (reading/writing ROOTfiles)
  TFile *inROOT;
  TFile *outROOT;

  //Input ROOTfile Name (to be read)
  TString data_InputFileName;
  TString data_InputReport;

  
  //Output ROOTfile Name
  TString data_OutputFileName;
  
  //ROOTfile to store combined hists from different runs
  TString data_OutputFileName_combined; 

  
  //Input parameter controls filenames
  TString main_controls_fname;
  TString input_CutFileName;
  TString input_HBinFileName;
  TString input_FileNamePattern;

  ///Output .txt filenames
  TString output_SummaryFileName;
  TString output_ReportFileName;

  //FileStreams objects to READ/WRITE to a .txt file
  ofstream out_file;
  ifstream in_file;

  // Get Counts of "total", "reals" and "random" events for saving to CaFe Report File"
  Double_t total_bins;

  Double_t W_total, W_rand, W_real, W_total_rate, W_real_rate;
  Double_t W_total_err, W_rand_err, W_real_err;
  
  Double_t Pm_total, Pm_rand, Pm_real, Pm_real_rate;
  Double_t Pm_total_err, Pm_rand_err, Pm_real_err;

  Double_t Em_total, Em_rand, Em_real;
  Double_t Em_total_err, Em_rand_err, Em_real_err;

  Double_t Em_nuc_total, Em_nuc_rand, Em_nuc_real;
  Double_t Em_nuc_total_err, Em_nuc_rand_err, Em_nuc_real_err;

  Double_t MM_total, MM_rand, MM_real;
  Double_t MM_total_err, MM_rand_err, MM_real_err;
  
  
  //------------DECLARE HISTOGRAM BINNING VARIABLES-------------
  Double_t nbins;

  //-------------------------
  //Particle Identification
  //-------------------------
  
  //Coincidence Time
  Double_t coin_nbins;
  Double_t coin_xmin;
  Double_t coin_xmax;

  //-HMS-
  Double_t hcer_nbins;	  
  Double_t hcer_xmin;	  
  Double_t hcer_xmax;	  
	                  
  Double_t hcal_nbins;	  
  Double_t hcal_xmin;	  
  Double_t hcal_xmax;	  
	                  
  Double_t hbeta_nbins;	  
  Double_t hbeta_xmin;	  
  Double_t hbeta_xmax;	  
  	                  
  //-SHMS-                
  Double_t pngcer_nbins;  
  Double_t pngcer_xmin;	  
  Double_t pngcer_xmax;	  
	                  
  Double_t phgcer_nbins;  
  Double_t phgcer_xmin;	  
  Double_t phgcer_xmax;	  
  	                  
  Double_t pcal_nbins;	  
  Double_t pcal_xmin;	  
  Double_t pcal_xmax;	  
	                  
  Double_t pbeta_nbins;	  
  Double_t pbeta_xmin;	  
  Double_t pbeta_xmax;	  
	                      

  //-----------------------------
  //Kinematics Histograms Bins
  //-----------------------------
  
  //Primary Kinematics
  Double_t the_nbins;
  Double_t the_xmin;
  Double_t the_xmax;

  Double_t W_nbins;
  Double_t W_xmin;
  Double_t W_xmax;

  Double_t W2_nbins;
  Double_t W2_xmin;
  Double_t W2_xmax;

  Double_t Q2_nbins;
  Double_t Q2_xmin;
  Double_t Q2_xmax;

  Double_t X_nbins;
  Double_t X_xmin;
  Double_t X_xmax;
  
  Double_t nu_nbins;
  Double_t nu_xmin;
  Double_t nu_xmax;

  Double_t q_nbins;
  Double_t q_xmin;
  Double_t q_xmax;

  Double_t qx_nbins;
  Double_t qx_xmin;
  Double_t qx_xmax;

  Double_t qy_nbins;
  Double_t qy_xmin;
  Double_t qy_xmax;

  Double_t qz_nbins;
  Double_t qz_xmin;
  Double_t qz_xmax;

  Double_t thq_nbins;
  Double_t thq_xmin;
  Double_t thq_xmax;

  Double_t phq_nbins;
  Double_t phq_xmin;
  Double_t phq_xmax;

  Double_t epsilon_nbins;
  Double_t epsilon_xmin;
  Double_t epsilon_xmax;

  //Secondary Kinematics

  Double_t Em_nbins;
  Double_t Em_xmin;
  Double_t Em_xmax;

  //nuclear missing energy (Mp + Mn - MA = nu - Tp - T_A-1 )
  Double_t Em_nuc_nbins;
  Double_t Em_nuc_xmin;
  Double_t Em_nuc_xmax;
  
  Double_t Pm_nbins;
  Double_t Pm_xmin;
  Double_t Pm_xmax;

  Double_t Pmx_lab_nbins;
  Double_t Pmx_lab_xmin;
  Double_t Pmx_lab_xmax;

  Double_t Pmy_lab_nbins;
  Double_t Pmy_lab_xmin;
  Double_t Pmy_lab_xmax;

  Double_t Pmz_lab_nbins;
  Double_t Pmz_lab_xmin;
  Double_t Pmz_lab_xmax;

  Double_t Pmx_q_nbins;
  Double_t Pmx_q_xmin;
  Double_t Pmx_q_xmax;

  Double_t Pmy_q_nbins;
  Double_t Pmy_q_xmin;
  Double_t Pmy_q_xmax;

  Double_t Pmz_q_nbins;
  Double_t Pmz_q_xmin;
  Double_t Pmz_q_xmax;

  Double_t Tx_nbins;
  Double_t Tx_xmin;
  Double_t Tx_xmax;

  Double_t Tr_nbins;
  Double_t Tr_xmin;
  Double_t Tr_xmax;
  
  Double_t MM_nbins;
  Double_t MM_xmin;
  Double_t MM_xmax;

  Double_t thxq_nbins;
  Double_t thxq_xmin;
  Double_t thxq_xmax;

  Double_t thrq_nbins;
  Double_t thrq_xmin;
  Double_t thrq_xmax;

  Double_t phxq_nbins;
  Double_t phxq_xmin;
  Double_t phxq_xmax;

  Double_t phrq_nbins;
  Double_t phrq_xmin;
  Double_t phrq_xmax;
  
  Double_t Tx_cm_nbins;
  Double_t Tx_cm_xmin;
  Double_t Tx_cm_xmax;

  Double_t Tr_cm_nbins;
  Double_t Tr_cm_xmin;
  Double_t Tr_cm_xmax;

  Double_t thxq_cm_nbins;
  Double_t thxq_cm_xmin;
  Double_t thxq_cm_xmax;

  Double_t thrq_cm_nbins;
  Double_t thrq_cm_xmin;
  Double_t thrq_cm_xmax;

  Double_t phxq_cm_nbins;
  Double_t phxq_cm_xmin;
  Double_t phxq_cm_xmax;

  Double_t phrq_cm_nbins;
  Double_t phrq_cm_xmin;
  Double_t phrq_cm_xmax;

  Double_t Ttot_cm_nbins;
  Double_t Ttot_cm_xmin;
  Double_t Ttot_cm_xmax;

  Double_t MandelS_nbins;
  Double_t MandelS_xmin;
  Double_t MandelS_xmax;

  Double_t MandelT_nbins;
  Double_t MandelT_xmin;
  Double_t MandelT_xmax;

  Double_t MandelU_nbins;
  Double_t MandelU_xmin;
  Double_t MandelU_xmax;

  //Kinematics Defined in HCANA (which are not in primary/secondary modules)
  Double_t kf_nbins;   //final electron momentum
  Double_t kf_xmin;
  Double_t kf_xmax;

  Double_t Pf_nbins;
  Double_t Pf_xmin;
  Double_t Pf_xmax;
  
  //Additional Kinematics
  Double_t thx_nbins;  //proton(hadron) angle
  Double_t thx_xmin;
  Double_t thx_xmax;

  Double_t MM2_nbins;  
  Double_t MM2_xmin;
  Double_t MM2_xmax;
  
  //----------------------------
  //Acceptance Histograms Bins
  //----------------------------


  //----Electron Arm Focal Plane-----
  Double_t exfp_nbins;
  Double_t exfp_xmin;
  Double_t exfp_xmax;

  Double_t expfp_nbins;
  Double_t expfp_xmin;
  Double_t expfp_xmax;
  
  Double_t eyfp_nbins;
  Double_t eyfp_xmin;
  Double_t eyfp_xmax;

  Double_t eypfp_nbins;
  Double_t eypfp_xmin;
  Double_t eypfp_xmax;

  //----Electron Arm Reconstructed-----
  Double_t eytar_nbins;
  Double_t eytar_xmin;
  Double_t eytar_xmax;

  Double_t eyptar_nbins;
  Double_t eyptar_xmin;
  Double_t eyptar_xmax;
  
  Double_t exptar_nbins;
  Double_t exptar_xmin ;
  Double_t exptar_xmax;
  
  Double_t edelta_nbins;
  Double_t edelta_xmin;  
  Double_t edelta_xmax;   

  
  //----Hadron Arm Focal Plane-----
  Double_t hxfp_nbins;
  Double_t hxfp_xmin;
  Double_t hxfp_xmax;

  Double_t hxpfp_nbins;
  Double_t hxpfp_xmin;
  Double_t hxpfp_xmax;
  
  Double_t hyfp_nbins;
  Double_t hyfp_xmin;
  Double_t hyfp_xmax;

  Double_t hypfp_nbins;
  Double_t hypfp_xmin;
  Double_t hypfp_xmax;

  //----Hadron Arm Reconstructed-----
  Double_t hytar_nbins;
  Double_t hytar_xmin;
  Double_t hytar_xmax;

  Double_t hyptar_nbins;
  Double_t hyptar_xmin;
  Double_t hyptar_xmax;
  
  Double_t hxptar_nbins;
  Double_t hxptar_xmin ;
  Double_t hxptar_xmax;
  
  Double_t hdelta_nbins;
  Double_t hdelta_xmin;  
  Double_t hdelta_xmax;   
  
  //----Target Quantities----
  //(Use same binning for hadron/electron reconstructed at target)
  Double_t tarx_nbins;
  Double_t tarx_xmin;
  Double_t tarx_xmax;

  Double_t tary_nbins;
  Double_t tary_xmin;
  Double_t tary_xmax;

  Double_t tarz_nbins;
  Double_t tarz_xmin;
  Double_t tarz_xmax;

  Double_t ztar_diff_nbins;
  Double_t ztar_diff_xmin;
  Double_t ztar_diff_xmax;
  
  //----Collimator Quantities----
  Double_t hXColl_nbins;
  Double_t hXColl_xmin;  
  Double_t hXColl_xmax;   
  
  Double_t hYColl_nbins;                                           
  Double_t hYColl_xmin;                                                                                                  
  Double_t hYColl_xmax;
  
  Double_t eXColl_nbins;
  Double_t eXColl_xmin;
  Double_t eXColl_xmax;
  
  Double_t eYColl_nbins;      
  Double_t eYColl_xmin;                                                                      
  Double_t eYColl_xmax;
  
  
  
  //-----------END SET DEFAULT HISTOGRAM BINNING----------


  //------------CREATE HISTOGRAMS--------------

  //Create Dummy histograms to store the ith and cumulative histograms (See CombineHistos() Method)
  //1D
  TH1F* h_total = 0;     //dummy histo to store hsitogram sum
  TH1F* h_i = 0;         //dummy histo to store ith histogram from list
  //2D
  TH2F* h2_total = 0;   //dummy histo to store hsitogram sum
  TH2F* h2_i = 0;       //dummy histo to store ith histogram from list

  //----------------------------------------------------------------
  // Detector Histograms (DATA ONLY): PID / TRACKING EFFICIENCY 
  //----------------------------------------------------------------
  
  //Coin. Time
  TH1F *H_ep_ctime_total;
  TH1F *H_ep_ctime;

  //HMS
  TH1F *H_hCerNpeSum;  
  TH1F *H_hCalEtotNorm;
  TH1F *H_hCalEtotTrkNorm;
  TH1F *H_hHodBetaNtrk;   
  TH1F *H_hHodBetaTrk;   
  
  //SHMS	       
  TH1F *H_pNGCerNpeSum;
  TH1F *H_pHGCerNpeSum;
  TH1F *H_pCalEtotNorm;
  TH1F *H_pCalEtotTrkNorm;
  TH1F *H_pHodBetaNtrk;   
  TH1F *H_pHodBetaTrk;   

  //-------Define 2D PID Histograms (correlations between pid detectors)-------
  //All possible combinations of correlations serves for general PID purposes

  //HMS 2D PID
  TH2F *H_hcal_vs_hcer;   //calorimeter vs. gas cherenkov
  
  //sHMS 2D PID
  TH2F *H_pcal_vs_phgcer;   //calorimeter vs. heavy gas cherenkov
  TH2F *H_pcal_vs_pngcer;   //calorimeter vs. noble gas cherenkov
  TH2F *H_pngcer_vs_phgcer;   //noble gas vs. heavy gas cherenkovs
  
  //DATA/SIMC Histograms (MUST BE THE EXACT SAME HSITOGRAM BINNING)

  //-----------------------------
  //   Kinematics Histograms
  //-----------------------------
  
  //Primary (electron) Kinematics
  TH1F *H_the;
  TH1F *H_W;
  TH1F *H_W2;
  TH1F *H_Q2;
  TH1F *H_xbj;
  TH1F *H_nu;
  TH1F *H_q;
  TH1F *H_qx;
  TH1F *H_qy;
  TH1F *H_qz;
  TH1F *H_thq;  
  TH1F *H_phq;   
  TH1F *H_epsilon;
  
  //Secondary (Hadron) Kinematics
  TH1F *H_Em;
  TH1F *H_Em_nuc;
  TH1F *H_Em_src;
  TH1F *H_Pm;
  TH1F *H_Pmx_lab;
  TH1F *H_Pmy_lab;
  TH1F *H_Pmz_lab;
  TH1F *H_Pmx_q;
  TH1F *H_Pmy_q;
  TH1F *H_Pmz_q;
  TH1F *H_Tx;
  TH1F *H_Tr;
  TH1F *H_MM;
  TH1F *H_thxq;
  TH1F *H_thrq;
  TH1F *H_phxq;
  TH1F *H_phrq;
  TH1F *H_Tx_cm;
  TH1F *H_Tr_cm;
  TH1F *H_thxq_cm;
  TH1F *H_thrq_cm;
  TH1F *H_phxq_cm;
  TH1F *H_phrq_cm;
  TH1F *H_Ttot_cm;
  TH1F *H_MandelS;
  TH1F *H_MandelT;
  TH1F *H_MandelU;

  //Kinematics Defined in HCANA (which are not in primary/secondary modules)
  TH1F *H_kf;
  TH1F *H_Pf;
  
  //Additional Kinematics (User-defined)
  TH1F *H_thx;    
  TH1F *H_MM2;     

  //----(Cosine, Sine) Histos of detecte AND recoil angles-----
  //(No need to define bins, as they are limited to -1 to 1). Will explicitly define bine in code
  //cth: cos(theta),  sth: sin(theta), cphi: cos(phi), sphi: sin(phi)
  TH1F *H_cth_xq;
  TH1F *H_cth_rq;
  TH1F *H_sth_xq;
  TH1F *H_sth_rq;
  TH1F *H_cphi_xq;
  TH1F *H_cphi_rq;
  TH1F *H_sphi_xq;
  TH1F *H_sphi_rq;

  TH1F *H_cth_xq_cm;
  TH1F *H_cth_rq_cm;
  TH1F *H_sth_xq_cm;
  TH1F *H_sth_rq_cm;
  TH1F *H_cphi_xq_cm;
  TH1F *H_cphi_rq_cm;
  TH1F *H_sphi_xq_cm;
  TH1F *H_sphi_rq_cm;

  // 2D Kinematics Histos
  TH2F *H_Em_nuc_vs_Pm;
  TH2F *H_Em_src_vs_Pm;
  
  //------------------------------------------------
  
  //-----------------------------
  //   Acceptance Histograms
  //-----------------------------
  
  //Electron Arm Focal Plane Quantities 
  TH1F *H_exfp;   
  TH1F *H_expfp;  
  TH1F *H_eyfp;   
  TH1F *H_eypfp;  

  //Electron Arm Reconstructed Quantities
  TH1F *H_eytar;  
  TH1F *H_eyptar; 
  TH1F *H_exptar; 
  TH1F *H_edelta; 

  //Hadron Arm Focal Plane Quantities 
  TH1F *H_hxfp;   
  TH1F *H_hxpfp;  
  TH1F *H_hyfp;   
  TH1F *H_hypfp;  
  //Hadron Arm Reconstructed Quantities
  TH1F *H_hytar;  
  TH1F *H_hyptar; 
  TH1F *H_hxptar; 
  TH1F *H_hdelta; 
  
  //Target Quantities (tarx, tary, tarz) in Hall Coord. System 
  TH1F *H_htar_x;  
  TH1F *H_htar_y;  
  TH1F *H_htar_z;  
  
  TH1F *H_etar_x;  
  TH1F *H_etar_y;  
  TH1F *H_etar_z;  

  //Additional Acceptance Quantities (User-defined)
  TH1F *H_ztar_diff;
  
  //Collimator Quantities (in spectrometer corrdinate system)
  TH1F *H_hXColl;  
  TH1F *H_hYColl;  
  TH1F *H_eXColl;  
  TH1F *H_eYColl;  
  
  //--2D Acceptance Histos--

  //Collimator Shape
  TH2F *H_hXColl_vs_hYColl;
  TH2F *H_eXColl_vs_eYColl;

  //Hour-Glass Shape
  TH2F *H_hxfp_vs_hyfp;
  TH2F *H_exfp_vs_eyfp;
  

  // ------- Selected Histograms for Random Coincidence Background Subtraction --------

  // NOTE: Nomenclature clarification 
  // *_rand -> random coincidence selection (sample selected outside the main coin. peak taken
  // to be representative of randoms underneath main coin. peak)  
  // *_rand_sub -> "true" coincidences after having subtracted the estimated randoms beneath the main peak

  //Coin. Time
  TH1F *H_ep_ctime_rand;
  TH1F *H_ep_ctime_rand_sub;

  // invariant mass
  TH1F *H_W_rand;
  TH1F *H_W_rand_sub;

  // 4-momentum transfer
  TH1F *H_Q2_rand;
  TH1F *H_Q2_rand_sub;

  // x-Bjorken
  TH1F *H_xbj_rand;
  TH1F *H_xbj_rand_sub;

  // energy transfer
  TH1F *H_nu_rand;
  TH1F *H_nu_rand_sub;

  // 3-momentum |q|
  TH1F *H_q_rand;
  TH1F *H_q_rand_sub;

  // missing energy 
  TH1F *H_Em_rand;
  TH1F *H_Em_rand_sub;

  // nuclear missing energy
  TH1F *H_Em_nuc_rand;
  TH1F *H_Em_nuc_rand_sub;

  // missing momentum
  TH1F *H_Pm_rand;
  TH1F *H_Pm_rand_sub;

  // missing mass
  TH1F *H_MM_rand;
  TH1F *H_MM_rand_sub;

  // in-plane angle between detected hadron and |q|
  TH1F *H_thxq_rand;
  TH1F *H_thxq_rand_sub;

  // in-plane angle between residual nucleus and |q|
  TH1F *H_thrq_rand;
  TH1F *H_thrq_rand_sub;

  
  //-----------END CREATE HISTOGRAMS-----------



  //----------DATA-RELATED VARIABLES-----------
  TTree *tree;
  Long64_t nentries;


  //Set-Up Tdc Counters for accepted triggers
  Double_t total_trig1_accp = 0;
  Double_t total_trig2_accp = 0;
  Double_t total_trig3_accp = 0;
  Double_t total_trig4_accp = 0;
  Double_t total_trig5_accp = 0;
  Double_t total_trig6_accp = 0;
  Double_t total_edtm_accp = 0;
  
  
  //Set-Up Tdc Counters for accepted triggers (passed bcm cut)
  Double_t total_trig1_accp_bcm_cut = 0;
  Double_t total_trig2_accp_bcm_cut = 0;
  Double_t total_trig3_accp_bcm_cut = 0;
  Double_t total_trig4_accp_bcm_cut = 0;
  Double_t total_trig5_accp_bcm_cut = 0;
  Double_t total_trig6_accp_bcm_cut = 0;
  Double_t total_edtm_accp_bcm_cut = 0;
  Double_t total_trig_accp_bcm_cut; //generic acc. trig

  //Tracking Efficiency Counter / Live Time (passed bcm cuts)
  //HMS
  Double_t h_did = 0;
  Double_t h_should = 0;
  Double_t hTrkEff;
  Double_t hTrkEff_err;
  
  //SHMS
  Double_t p_did = 0;
  Double_t p_should = 0;
  Double_t pTrkEff;
  Double_t pTrkEff_err;

  //Computer Live Time 
  Double_t cpuLT_trig;       //generic computer live time
  Double_t cpuLT_trig_err_Bi;  //generic cpu live time error (using binomial statistics)
  Double_t cpuLT_trig_err_Bay;  //generic cpu live time error (using bayesian statistics)
   
  Double_t cpuLT_trig1, cpuLT_trig1_err_Bi, cpuLT_trig1_err_Bay;
  Double_t cpuLT_trig2, cpuLT_trig2_err_Bi, cpuLT_trig2_err_Bay;
  Double_t cpuLT_trig3, cpuLT_trig3_err_Bi, cpuLT_trig3_err_Bay;
  Double_t cpuLT_trig4, cpuLT_trig4_err_Bi, cpuLT_trig4_err_Bay;
  Double_t cpuLT_trig5, cpuLT_trig5_err_Bi, cpuLT_trig5_err_Bay;
  Double_t cpuLT_trig6, cpuLT_trig6_err_Bi, cpuLT_trig6_err_Bay;

  //Total Live Time (EDTM)
  Double_t tLT_corr_factor;

  Double_t tLT_trig; //generic total live time 
  Double_t tLT_trig_err_Bi;  //generic total live time error (using binomial statistics)
  Double_t tLT_trig_err_Bay;  //generic total live time error (using bayesian statistics)

  Double_t tLT_trig1, tLT_trig1_err_Bi, tLT_trig1_err_Bay;
  Double_t tLT_trig2, tLT_trig2_err_Bi, tLT_trig2_err_Bay;
  Double_t tLT_trig3, tLT_trig3_err_Bi, tLT_trig3_err_Bay;
  Double_t tLT_trig4, tLT_trig4_err_Bi, tLT_trig4_err_Bay;
  Double_t tLT_trig5, tLT_trig5_err_Bi, tLT_trig5_err_Bay;
  Double_t tLT_trig6, tLT_trig6_err_Bi, tLT_trig6_err_Bay;

  
  //----------DEFINE CUT VARIABLES USED IN LIVE TIME / TRACKING EFFICIENCY DEFINITIONS (DATA ONLY)---------

  //edtm cuts
  Bool_t c_noedtm;
  Bool_t c_edtm;

  //generic cpu live time cuts (HMS singles -> htrig, SHMS singles -> ptrig, coin -> ptrig)
  Bool_t c_trig1;    
  Bool_t c_trig2;   
  Bool_t c_trig3;   
  Bool_t c_trig4;   
  Bool_t c_trig5;   
  Bool_t c_trig6;

  //Pre-Scale factor for each pre-trigger (used in computer/total live time calculation, to account for pre-scaled triggers)
  Float_t Ps_factor = 1;  //generic pre-scale factor
  Float_t Ps1_factor = 1;  //default is 1
  Float_t Ps2_factor = 1;
  Float_t Ps3_factor = 1;
  Float_t Ps4_factor = 1;
  Float_t Ps5_factor = 1;
  Float_t Ps6_factor = 1; 
  
  //PID cuts (for tracking efficiency. If used for data anylsis pid, it will be the exact same as pid used for tracking efficiency, so it will not be biased)

  //HMS
  Bool_t hdc_ntrk_cut_flag;    Bool_t c_hdc_ntrk;       Double_t c_hdc_ntrk_min;
  Bool_t hScinGood_cut_flag;   Bool_t c_hScinGood;  // on or off (no limit required)
  Bool_t hcer_cut_flag;        Bool_t c_hcer_NPE_Sum;   Double_t c_hnpeSum_min;    Double_t c_hnpeSum_max;    
  Bool_t hetotnorm_cut_flag;   Bool_t c_hetotnorm;      Double_t c_hetotnorm_min;  Double_t c_hetotnorm_max;
  Bool_t hBeta_notrk_cut_flag; Bool_t c_hBeta_notrk;    Double_t c_hBetaNtrk_min;  Double_t c_hBetaNtrk_max;   
  
  //SHMS
  Bool_t pdc_ntrk_cut_flag;    Bool_t c_pdc_ntrk;         Double_t c_pdc_ntrk_min;
  Bool_t pScinGood_cut_flag;   Bool_t c_pScinGood;   // on or off (no limit required)

  Bool_t pngcer_cut_flag;      Bool_t c_pngcer_NPE_Sum;   Double_t c_pngcer_npeSum_min;  Double_t c_pngcer_npeSum_max;    
  Bool_t phgcer_cut_flag;      Bool_t c_phgcer_NPE_Sum;   Double_t c_phgcer_npeSum_min;  Double_t c_phgcer_npeSum_max;   
  Bool_t petotnorm_cut_flag;   Bool_t c_petotnorm;        Double_t c_petotnorm_min;      Double_t c_petotnorm_max;
  Bool_t pBeta_notrk_cut_flag; Bool_t c_pBeta_notrk;      Double_t c_pBetaNtrk_min;      Double_t c_pBetaNtrk_max;

  
  //e- tracking efficiency Boolean
  Bool_t good_hms_should;
  Bool_t good_hms_did;
  
  //h tracking efficiency Boolean
  Bool_t good_shms_should;
  Bool_t good_shms_did;
  
  
  //-----------DEFINE DATA/SIMC CUTS (CUTS MUST BE EXACT SAME. Which is why only a variable is used for both)------
  //(See set_basic_cuts.inp file to modify the cuts),  the 'c_' denotes it is a cut

  //------STANDARD PID Cuts ON DATA (THESE ARE SLIGHTLY DIFFERENT THAN IN TRK EFF. DEFINITION -- SHOULD THEY BE THE SAME??)

  //Coincidence time cut 
  Bool_t ePctime_cut_flag;
  Bool_t eP_ctime_cut;
  Bool_t eP_ctime_cut_rand; 
  Bool_t eP_ctime_cut_rand_L;  // boolean for selecting random coincidences left (outside main coin. peak)
  Bool_t eP_ctime_cut_rand_R;  // boolean for selecting random coincidences right (outside main coin. peak)
  //Double_t ePctime_cut_thrs; // coin time threshold cut: coin_time_peak +/- cpid_ePctime_thrs

  // main coincidence time peak min/max window cut
  Double_t ePctime_cut_min;
  Double_t ePctime_cut_max;
  Double_t dt_coin_peak;  // main coin. time window width
  
  // accidentals to the right of main coin. peak
  Double_t ePctime_cut_min_R;
  Double_t ePctime_cut_max_R;
  Double_t dt_acc_R;   // window width of accidentals (right of main coin. peak)
  
  // accidentals to the left of main coin. peak
  Double_t ePctime_cut_min_L;
  Double_t ePctime_cut_max_L;
  Double_t dt_acc_L;   // window width of accidentals (left of main coin. peak)
  
  //coin time integer multiple of eP_ctime_thrs cut used to select randoms (must be: >=2: i,e, 2, 3, 4, . . .) 
  Float_t eP_mult;           

  //Scale factor variables (for random coincidence scaling / subtraction)
  Float_t P_scale_factor;

  //SHMS Calorimeter EtotTrackNorm (e- selection)
  Bool_t petot_trkNorm_pidCut_flag;
  Bool_t cpid_petot_trkNorm;
  Double_t cpid_petot_trkNorm_min;
  Double_t cpid_petot_trkNorm_max;

  //SHMS Noble Gas Cherenkov
  Bool_t pngcer_pidCut_flag;
  Bool_t cpid_pngcer_NPE_Sum;
  Double_t cpid_pngcer_npeSum_min;
  Double_t cpid_pngcer_npeSum_max;    

  //SHMS Heavy Gas Cherenkov (Pi /K separation)
  Bool_t phgcer_pidCut_flag;
  Bool_t cpid_phgcer_NPE_Sum;
  Double_t cpid_phgcer_npeSum_min;
  Double_t cpid_phgcer_npeSum_max;

  //HMS Calorimeter EtotTrackNorm (e- selection)
  Bool_t hetot_trkNorm_pidCut_flag;
  Bool_t cpid_hetot_trkNorm;
  Double_t cpid_hetot_trkNorm_min;
  Double_t cpid_hetot_trkNorm_max;

  //HMS Gas Cherenkov
  Bool_t hcer_pidCut_flag;
  Bool_t cpid_hcer_NPE_Sum;
  Double_t cpid_hcer_npeSum_min;
  Double_t cpid_hcer_npeSum_max;    
  
  //---------Kinematics Cuts----------

  // H(e,e'p) Kinematic Cuts
  
  //4-Momentum Transfers
  Bool_t Q2_heep_cut_flag;
  Bool_t c_heep_Q2;
  Double_t c_heep_Q2_min;
  Double_t c_heep_Q2_max;

  //4-Momentum Transfers
  Bool_t xbj_heep_cut_flag;
  Bool_t c_heep_xbj;
  Double_t c_heep_xbj_min;
  Double_t c_heep_xbj_max;
  
  //Missing Energy
  Bool_t Em_heep_cut_flag;
  Bool_t c_heep_Em;
  Double_t c_heep_Em_min;
  Double_t c_heep_Em_max;

  //Invariant Mass, W
  Bool_t W_heep_cut_flag;
  Bool_t c_heep_W;
  Double_t c_heep_W_min;
  Double_t c_heep_W_max;

  //Missing Mass, MM
  //Protons
  Bool_t MM_heep_cut_flag;
  Bool_t c_heep_MM;
  Double_t c_heep_MM_min;
  Double_t c_heep_MM_max;

  // CaFe A(e,e'p) Mean-Field (MF) Kinematic Cuts -----
  Bool_t   Q2_MF_cut_flag;
  Bool_t   c_MF_Q2;
  Double_t c_MF_Q2_min;
  Double_t c_MF_Q2_max;

  Bool_t   Pm_MF_cut_flag;
  Bool_t   c_MF_Pm;
  Double_t c_MF_Pm_min;
  Double_t c_MF_Pm_max;

  // missing energy cut (only for deuteron, since Em ~ 2.2 MeV for deuteron)
  Bool_t   Em_d2MF_cut_flag;
  Bool_t   c_d2MF_Em;
  Double_t c_d2MF_Em_min;
  Double_t c_d2MF_Em_max;
  
  // missing energy cut (only for A>2 nuclei)
  Bool_t   Em_MF_cut_flag;
  Bool_t   c_MF_Em;
  Double_t c_MF_Em_min;
  Double_t c_MF_Em_max;
  
  // CaFe A(e,e'p) Short-Range Correlations (SRC) Kinematic Cuts -----
  Bool_t   Q2_SRC_cut_flag;
  Bool_t   c_SRC_Q2;
  Double_t c_SRC_Q2_min;
  Double_t c_SRC_Q2_max;

  Bool_t   Pm_SRC_cut_flag;
  Bool_t   c_SRC_Pm;
  Double_t c_SRC_Pm_min;
  Double_t c_SRC_Pm_max;

  Bool_t   Xbj_SRC_cut_flag;
  Bool_t   c_SRC_Xbj;
  Double_t c_SRC_Xbj_min;
  Double_t c_SRC_Xbj_max;
  
  Bool_t   thrq_SRC_cut_flag;
  Bool_t   c_SRC_thrq;
  Double_t c_SRC_thrq_min;
  Double_t c_SRC_thrq_max;

  Bool_t   Em_d2SRC_cut_flag;
  Bool_t   c_d2SRC_Em;
  Double_t c_d2SRC_Em_min;
  Double_t c_d2SRC_Em_max;

  // Missing energy cut on A>2 nuclei for SRC kinematics (this is a dynamic cut, which varies with Pm)
  // is itslef a square-root function of Pm, Em_src (Pm) = nu - Tp - (sqrt(Pm*Pm + MN*MN) - MN)
  Bool_t   Em_SRC_cut_flag;
  Bool_t   c_SRC_Em;
  
  //----------Acceptance Cuts------------  

  //-Hadron Arm-

  // Momentum Acceptance, Delta [%] |
  Bool_t hdelta_cut_flag;
  Bool_t c_hdelta;
  Double_t c_hdelta_min;
  Double_t c_hdelta_max;

  // Angular Acceptance
  Bool_t hxptar_cut_flag;
  Bool_t c_hxptar;
  Double_t c_hxptar_min;
  Double_t c_hxptar_max;

  Bool_t hyptar_cut_flag;
  Bool_t c_hyptar;
  Double_t c_hyptar_min;
  Double_t c_hyptar_max;
  
  //-Electron Arm-

  // Momentum Acceptance, Delta [%] 
  Bool_t edelta_cut_flag;
  Bool_t c_edelta;
  Double_t c_edelta_min;
  Double_t c_edelta_max;

  // Angular Acceptance
  Bool_t exptar_cut_flag;
  Bool_t c_exptar;
  Double_t c_exptar_min;
  Double_t c_exptar_max;

  Bool_t eyptar_cut_flag;
  Bool_t c_eyptar;
  Double_t c_eyptar_min;
  Double_t c_eyptar_max;
  
  // Z-Reaction Vertex Difference Cut
  Bool_t ztarDiff_cut_flag;
  Bool_t c_ztarDiff;
  Double_t c_ztarDiff_min;
  Double_t c_ztarDiff_max;
  
  //COMBINE ALL CUTS
  Bool_t c_baseCuts;  //base cuts

  Bool_t c_accpCuts_hms;  //acceptance cuts
  Bool_t c_accpCuts_shms;  //acceptance cuts
  Bool_t c_accpCuts;  //acceptance cuts

  Bool_t c_pidCuts_hms;   //particle id cuts
  Bool_t c_pidCuts_shms;   //particle id cuts
  Bool_t c_pidCuts;   //particle id cuts
  
  Bool_t c_kinHeepSing_Cuts;     //kinematics cuts (Heep e- Singles Cuts)
  Bool_t c_kinHeepCoin_Cuts;     //kinematics cuts (Heep Coin Cuts)
  Bool_t c_kinMF_Cuts;  //kinematics cuts (CaFe MF Cuts)
  Bool_t c_kinSRC_Cuts; //kinematics cuts (CaFe SRC Cuts)
  
  //------------------END DATA-RELATED VARIABLES DEFINED CUTS-------------------
  
  //-------------TTREE LEAF VARIABLE NAMES (DATA or SIMC)--------------

  //Trigger Detector / Global Variables
  Double_t gevtyp;  //global event type
  Double_t gevnum;
  Double_t TRIG1_tdcTimeRaw;
  Double_t TRIG2_tdcTimeRaw;
  Double_t TRIG3_tdcTimeRaw;
  Double_t TRIG4_tdcTimeRaw;
  Double_t TRIG5_tdcTimeRaw;
  Double_t TRIG6_tdcTimeRaw;
  Double_t EDTM_tdcTimeRaw;

  //--------------------------------
  //-------PID Leaf Variables-------
  //--------------------------------
  //Coincidence Time (ONLY DATA)
  Double_t epCoinTime;
  Double_t eKCoinTime;
  Double_t ePiCoinTime;

  //HMS DETECTORS
  Double_t hcer_npesum;
  Double_t hcal_etotnorm;
  Double_t hcal_etottracknorm;
  Double_t hhod_beta_ntrk;
  Double_t hhod_beta;
  Double_t hhod_gtr_beta;
  Double_t hhod_GoodScinHit;
  Double_t hdc_ntrack;    
 
  //SHMS DETECTORS
  Double_t phgcer_npesum;
  Double_t pngcer_npesum;
  Double_t pcal_etotnorm;
  Double_t pcal_etottracknorm;
  Double_t phod_beta_ntrk;
  Double_t phod_beta;
  Double_t phod_gtr_beta;
  Double_t phod_GoodScinHit;
  Double_t pdc_ntrack;    

  //-----------------------------------
  //-----Kinematics Leaf Variables-----
  //-----------------------------------
  //Primary Kinematics (electron kinematics)
  Double_t th_e;                 //electron arm particle central angle relative to +z (hall coord. system) [rad]
  Double_t W;                    //Invariant Mass W assuming the proron mass [GeV]
  Double_t W2;                   //Invariant Mass W2 (squared) assuming the proron mass [GeV^2]  
  Double_t Q2;                   //Four-momentum trasfer [(GeV/c)^2]
  Double_t X;                    //B-jorken X  scaling variable
  Double_t nu;                   //Energy Transfer [GeV]
  Double_t q;                    //magnitude of the 3-vector q [GeV]
  Double_t qx;                   //x-comp of the 3-vector q (lab) 
  Double_t qy;                   //y-comp of the 3-vector q (lab)
  Double_t qz;                   //z-comp of the 3-vector q (lab)
  Double_t th_q;                 //in-plane angle between q and +z (hall coord. system, along beam direction) [rad]
  Double_t ph_q;                 //out-of-plane angle between q and +z (hall coord. system) [rad]
  Double_t epsilon;              //virtual photon polarization factor
  
  //Secondary Kinematics (USED BY DATA AND SIMC)
  //(missing and recoil are used inter-changeably)
  Double_t Em;                    //Standard Missing Energy for H(e,e'p) [GeV]
  Double_t Em_nuc;                //Nuclear definition of Missing Energy (Used for D(e,e'p): B.E. of deuteron) [GeV]

  Double_t Em_src;                // Em_src = nu - Tp - T_"n", where T_"n" = E_nucleon - M_nucleon = sqrt(M_nucleon^2 + Pm^2) - M_nucleon,
                                  // T_"n" is assumed to be kin. energy of second SRC nucleon. So basically, missing energy of recoil nucleon in SRC pair

  Double_t Pm;                    // Missing Momentum (should be zero for H(e,e'p) [GeV/c]
  Double_t Pmx_lab;               //X-Component of Missing Momentum (in Lab(or Hall) frame. +X: beam left, +Y: up, +Z: downstream beam) 
  Double_t Pmy_lab;
  Double_t Pmz_lab;
  Double_t Pmx_q;                 //X-Component of Missing Momentum (in frame where +z_lab is rotated to +z_q. Pmz_q is along +z(parallel to q))
  Double_t Pmy_q;                 //Y-component of Missing Momentum w.r.to q [GeV/c]
  Double_t Pmz_q;                 //Z-Component of Missing Momentum along +q [GeV/c]
  Double_t Tx;                    //Kinetic Energy of detected particle (proton, kaon, pion, . . .) [GeV]
  Double_t Tr;                    //Kinetic Energy of recoil system [GeV]
  Double_t MM;                    //Invariant ('Missing Mass') of recoil system [GeV/c^2]
  Double_t th_xq;                 //In-plane angle between detected particle and q [rad]  
  Double_t th_rq;                 //In-plane angle between the recoil system and q [rad]  
  Double_t ph_xq;                 //Out-of-plane angle between detected particle and q [rad]   
  Double_t ph_rq;                 //Out-of-plane anfle between recoil system and q [rad]
  Double_t xangle;                //Angle of detected particle with scattered electron (Used to determine hadron angle) [rad]
  Double_t Tx_cm;                 //Kinetic Energy of detected particle in CM frame (proton, kaon, pion, . . .) [GeV]
  Double_t Tr_cm;                 //Kinetic Energy of recoil system in CM frame [GeV]
  Double_t th_xq_cm;              //In-plane angle between detected particle and q in CM frame [rad]  
  Double_t th_rq_cm;              //In-plane angle between the recoil system and q in CM frame [rad]  
  Double_t ph_xq_cm;              //Out-of-plane angle between detected particle and q in CM frame [rad]   
  Double_t ph_rq_cm;              //Out-of-plane anfle between recoil system and q in CM frame [rad]
  Double_t Ttot_cm;               //Total kinetic energy in CM frame [GeV]
  Double_t MandelS;               //Mandelstam s for secondary vertex [GeV^2]
  Double_t MandelT;               //Mandelstam t for secondary vertex [GeV^2]
  Double_t MandelU;               //Mandelstam u for secondary vertex [GeV^2]

  //Kinematics Defined in HCANA (which are not in primary/secondary modules)
  Double_t kf;       // final electron arm momentum [GeV/c]
  Double_t Pf;       // final proton momentum [GeV/c]
  
  //Additional Kinematics (User-defined)
  Double_t th_x;                   //hadron arm particle central angle
  Double_t MM2;                   //Missing Mass Squared

  //------------------------------------
  //-----Acceptance Leaf Variables------
  //------------------------------------
  
  //Electron Arm Focal Plane Quantities 
  Double_t e_xfp;   // [cm]
  Double_t e_xpfp;  // [rad]
  Double_t e_yfp;   // [cm]
  Double_t e_ypfp;  // [rad]

  //Electron Arm Reconstructed Quantities
  Double_t e_ytar;   // [cm]
  Double_t e_yptar;  // [rad]
  Double_t e_xptar;  // [rad]
  Double_t e_delta;  // percent

  //Hadron Arm Focal Plane / Reconstructed Quantities (USED BY DATA AND SIMC)
  Double_t h_xfp;   // [cm]
  Double_t h_xpfp;  // [rad]
  Double_t h_yfp;   // [cm]
  Double_t h_ypfp;  // [rad]
  
  Double_t h_ytar;   // [cm]
  Double_t h_yptar;  // [rad]
  Double_t h_xptar;  // [rad]
  Double_t h_delta;  // percent
  
  //Target Quantities (tarx, tary, tarz) in Hall Coord. System (USED BY DATA AND SIMC)
  Double_t htar_x;  //[cm]
  Double_t htar_y;  //[cm]
  Double_t htar_z;  //[cm]
  
  Double_t etar_x;  //[cm]
  Double_t etar_y;  //[cm]
  Double_t etar_z;  //[cm]

  //Collimator Quantities (in spectrometer corrdinate system)
  Double_t hXColl;  //[cm]
  Double_t hYColl;  //[cm]
  Double_t eXColl;  //[cm]
  Double_t eYColl;  //[cm]
  
  //Additional Acceptance Quantities (User-defined)
  Double_t ztar_diff; //[cm]


  
  //----------END TTREE LEAF VARIABLE NAMES (DATA or SIMC)--------------


  //----------SCALER-RELATED VARIABLES----------
  TTree *scaler_tree;
  Long64_t scal_entries;
  Int_t *evt_flag_bcm;   //flag (0 or 1) to determine whether the scaler read passed the cut
  Int_t *scal_evt_num;   //store data event number associated with scaler read 
  Double_t set_current;  //Set current for each run. For now, take the current ->  maximum bin content of bcm current histogram 
  Int_t scal_read = 0;   //scaler read counter (actually used inside data loop to count scaler reads)

  //Define Counter Quantities To Store Previous Reads
  Double_t prev_time = 0.;
  Double_t prev_charge_bcm = 0.;
  Double_t prev_s1x_scaler = 0;
  Double_t prev_trig1_scaler = 0;
  Double_t prev_trig2_scaler = 0;
  Double_t prev_trig3_scaler = 0;
  Double_t prev_trig4_scaler = 0;
  Double_t prev_trig5_scaler = 0;
  Double_t prev_trig6_scaler = 0;
  Double_t prev_edtm_scaler = 0;

  //Define Counter Quantities To Store Accumulated Reads
  Double_t total_time = 0.;
  Double_t total_charge_bcm = 0.;
  Double_t total_s1x_scaler = 0;
  Double_t total_trig1_scaler = 0;
  Double_t total_trig2_scaler = 0;
  Double_t total_trig3_scaler = 0;
  Double_t total_trig4_scaler = 0;
  Double_t total_trig5_scaler = 0;
  Double_t total_trig6_scaler = 0;
  Double_t total_edtm_scaler = 0;

  //Store Accumulated Reads if they passed BCM Current Cut
  Double_t total_time_bcm_cut = 0.;
  Double_t total_charge_bcm_cut = 0.;
  Double_t total_s1x_scaler_bcm_cut = 0;
  Double_t total_trig1_scaler_bcm_cut = 0;
  Double_t total_trig2_scaler_bcm_cut = 0;
  Double_t total_trig3_scaler_bcm_cut = 0;
  Double_t total_trig4_scaler_bcm_cut = 0;
  Double_t total_trig5_scaler_bcm_cut = 0;
  Double_t total_trig6_scaler_bcm_cut = 0;
  Double_t total_edtm_scaler_bcm_cut = 0;
  Double_t total_trig_scaler_bcm_cut; //generic trig scaler count

  //Store Scaler Rates if current cut passed
  Double_t S1XscalerRate_bcm_cut;
  Double_t TRIG1scalerRate_bcm_cut;
  Double_t TRIG2scalerRate_bcm_cut;
  Double_t TRIG3scalerRate_bcm_cut;
  Double_t TRIG4scalerRate_bcm_cut;
  Double_t TRIG5scalerRate_bcm_cut;
  Double_t TRIG6scalerRate_bcm_cut;
  Double_t EDTMscalerRate_bcm_cut;

  //Store Accepted Rates if current cut passed
  Double_t TRIG1accpRate_bcm_cut;
  Double_t TRIG2accpRate_bcm_cut;
  Double_t TRIG3accpRate_bcm_cut;
  Double_t TRIG4accpRate_bcm_cut;
  Double_t TRIG5accpRate_bcm_cut;
  Double_t TRIG6accpRate_bcm_cut;
  Double_t EDTMaccpRate_bcm_cut;
  
  Double_t trig_rate; //generic trigger rate
  
  //Store Average BCM Current
  Double_t  avg_current_bcm_cut;
  
  //---------END SCALER-RELATED VARIABLES----------


  //--------SCALER TTREE VARIABLE NAMES (DATA)---------

  Double_t Scal_evNum;
  Double_t Scal_BCM_charge;
  Double_t Scal_BCM_current;
  Double_t Scal_time;
  Double_t S1X_scaler;
  Double_t TRIG1_scaler;   
  Double_t TRIG2_scaler;   
  Double_t TRIG3_scaler;   
  Double_t TRIG4_scaler;   
  Double_t TRIG5_scaler;   
  Double_t TRIG6_scaler;   
  Double_t EDTM_scaler;
  
  //-------END SCALER TTREE VARIABLE NAMES (DATA)------

  //-----------------VARIABLES RELATED TO FullWeight APPLIED TO DATA-YIELD------------------

  Double_t FullWeight = 1;  //default
  
  Double_t hadAbs_corr;           //correct for lost coincidences due to the hadron in HMS (or SHMS) NOT making it to the hodoscopes to form trigger
  Double_t hadAbs_corr_err;       //uncertainty in hadron absorption correction factor
  Double_t tgtBoil_corr;          //correct for lost coincidences due to localized target boiling at high currents
  Double_t tgtBoil_corr_err;      //uncertainty in target boiling correction factor
  
  //Slope of Norm. Yield vs. avg beam current (Determined from target boiling studies)
  Double_t tgt_slope, tgt_slope_err;  //generic target slope
  Double_t LH2_slope;     
  Double_t LD2_slope;
  Double_t LH2_slope_err;     
  Double_t LD2_slope_err;
  
  //beam current relative and absolute errors (necessary in error propagation calculation of tgt boil corr. factor)
  Double_t dIb_Ib;   
  Double_t avg_current_bcm_cut_err;


  //-----------------------------------------------------------------------------------------


  //------VARIABLES USED TO WRITE HISTOGRAMS TO ROOT FILE-------

  //Create Categorical TLists to store histograms based on caterogy
  TList *pid_HList;    //store detector histograms (i.e., coin_time, H_hcerNpeSum, H_eCalEtotNorm, . . .)
  TList *kin_HList;    //store kinematical histograms (i.e., Q2, W, th_e, . . .)
  TList *accp_HList;   //store spectrometer accpetance histograms (focal plane, reconstructed, ztar_diff)

  TList * rand_HList;  // store random coin. background of selected histograms
  TList * randSub_HList;  // store random-subtracted variables of selected histograms
  
  //---------------------------------------------




  //-------------------------
  // HELICITY TREE VARIABLES
  //-------------------------

  //These variables are meant to be used by the derived helicity class 
  //We will have a separate EventLoop for the helicity analysis
  //have NOT figure out yet how to make additions to a base class in a dervied class and ONLY call the
  //derived method once.
  
  //Data TTree Leaf Variables
  Double_t hel_cycle;      //Helicity Cycle
  Double_t hel;           //actual helicity for event
  Double_t hel_pred;       //predicted reported helicity for event
  Double_t hel_rep;        //reported helicity for event
  Double_t hel_mps;        //In MPS blanking period (helicity unknown or undetermined)
  Double_t hel_nqrt;       //position of helicity in quartet
  Double_t hel_pcheck;     //Period check
  Double_t hel_qrt;        //Last cycle of quartet
  
  
  
};

#endif 
