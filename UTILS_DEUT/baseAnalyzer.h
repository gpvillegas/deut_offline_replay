/*
Author: Carlos Yero
email: cyero002@fiu.edu, cyero@jlab.org
Date Created: August 22, 2020
*/

#ifndef BASE_ANALYZER_H
#define BASE_ANALYZER_H

#include "../UTILS/parse_utils.h" //useful C++ string parsing utilities
#include "../UTILS/hist_utils.h" //useful C++ histogram bin extraction utility
#include "../UTILS/project2d.h" //useful C++ 2d histogram projection utility

#include <string>

class baseAnalyzer
{

  
public:
  
  //Constructor / Destructor
  baseAnalyzer( int irun=-1, int ievt=-1, string mode="", string earm="", string ana_type="", string ana_cuts="", Bool_t hel_flag=0, string bcm_name="", double thrs=-1, string trig_single="", string trig_coin="", Bool_t create_skim=0); //initialize member variables

  //2nd constructor (overload constructor, i.e., different arguments)
  baseAnalyzer(string earm="", string ana_type="", string ana_cuts="");
  
  ~baseAnalyzer();
  
  //MAIN ANALYSIS FUNCTIONS
  void run_online_data_analysis();
  void run_offline_data_analysis();
  void run_simc_analysis();
  void run_scalers(); // mainly for generating cafe output file (for bcm calib runs)
  
  //Function prototypes
  void Init(); 
  void ReadInputFile();
  void ReadReport();
  void SetHistBins();
  void CreateHist();
  void ReadScalerTree();  
  void ScalerEventLoop(); //bcm current cut threshold in uA units
  void ReadTree();
  void CreateSkimTree();
  void CreateSinglesSkimTree();
  void EventLoop();
  void CalcEff();
  void ApplyWeight();
  void ScaleSIMC(TString target="");
  void WriteHist();
  void WriteReport();
  void WriteOnlineReport();
  void WriteOfflineReport();
  void WriteReportSummary();
  void CombineHistos();
  void TrackOnlineStats();
  
  //void CalcRadCorr(); 
  //void ApplyRadCorr();
  //void ChargeNorm(); 
  void RandSub(); //Apply subtraction of random coincidence background
  //void GetAsymmetry();
  
  // Helper Functions
  void GetPeak();
  void CollimatorStudy();
  void MakePlots();
  Double_t GetLuminosity(TString user_input="");

    
  //hms/shms dc calibration quality monitoring constanta
  static const Int_t dc_PLANES = 12;
  const string hdc_pl_names[dc_PLANES] = {"1u1", "1u2", "1x1", "1x2", "1v2", "1v1", "2v1", "2v2", "2x2", "2x1", "2u2", "2u1"};
  const string pdc_pl_names[dc_PLANES] = {"1u1", "1u2", "1x1", "1x2", "1v1", "1v2", "2v2", "2v1", "2x2", "2x1", "2u2", "2u1"};

protected:

  //Set Constants
  const Double_t pi = TMath::Pi(); 
  const Double_t dtr = pi / 180.;
  const Double_t amu2GeV = 0.93149432;
  const Double_t cm2topb = 1./(1e-36);   // 1 pb / (1e-36 cm^2)
  const Double_t pb_to_invGeV2 = 1. / (2.56819*1e-9); // 1 pb / (2.56819*1e-9 GeV^-2)
  const Double_t cm2_to_invGeV2 = 1. / (3.8937929*1e-28);    // (1e-36 cm^2/ 1 pb ) * 1 pb / (2.56819*1e-9 GeV^-2)
  const Double_t NA = 6.022*1e23;  // Avogadro's number ( # atoms / mol), 1 g/mol = 1amu
  const Double_t elementary_charge = 1.60217663*1e-19; // Coulombs 
  
  // target mass (amu),  mass number (# of nucleons)
  Int_t mass_number_A = 0;
  Double_t MH_amu     = 1.00794       , A_H   = 1;
  Double_t MD_amu     = 2.01410177812 , A_D   = 2;
  Double_t MBe9_amu   = 9.012182      , A_Be9 = 9;
  Double_t MB10_amu   = 10.0129370    , A_B10 = 10;   // target is actually 10B4C (Boron-Carbide)
  Double_t MB11_amu   = 11.009306     , A_B11 = 11;  // target is actually 11B4C (Boron-Carbide)
  Double_t MC12_amu   = 12.0107       , A_C12 = 12;
  Double_t MAl27_amu  = 26.98153      , A_Al27 = 27;
  Double_t MCa40_amu  = 39.962590863  , A_Ca40 = 40;
  Double_t MCa48_amu  = 47.95252276   , A_Ca48 = 48;
  Double_t MFe54_amu  = 53.9396147    , A_Fe54 = 54;
  Double_t MTi48_amu  = 47.9479463    , A_Ti48 = 48;
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
  
  //target length (thickness) (cm) obtained from D. Meekins table on Hall C 2022 targets
  //https://docs.google.com/spreadsheets/d/1GoHMHbCv3v6CbVybqTQws-4VhQeSifFP/edit#gid=140150071

  // # of nucleons
  Int_t N;  // neutrons
  Int_t Z;  // protons
  Int_t A;  // A = N+Z nucleons
  Double_t T = 0.0 ; //transparency placeholder
  
  //target density (g/cm^3)
  Double_t tgt_density = 0.0;  // generic variable to hold target density
  Double_t rho_H    = 0.07231;   
  Double_t rho_D    = 0.167;   
  Double_t rho_Be9  = 1.848; //
  Double_t rho_B10  = 2.352; // 
  Double_t rho_B11  = 2.434; // 
  Double_t rho_C12  = 1.8;  //
  Double_t rho_Al27 = 2.699; //
  Double_t rho_Ca40 = 1.55; //
  Double_t rho_Ca48 = 1.86; //
  Double_t rho_Fe54 = 7.87;
  Double_t rho_Ti48 = 4.5;

  //target thickness (or length) (cm)
  Double_t tgt_thickness = 0.0;  // generic variable to hold target thickness
  Double_t thick_H    = 10.;   
  Double_t thick_D    = 10.;
  Double_t thick_Be9  = 0.5335; //
  Double_t thick_B10  = 0.245; //
  Double_t thick_B11  = 0.26; //
  Double_t thick_C12  = 0.3188; //
  Double_t thick_Al27 = 0.0889; // 0.0889 (uptream), 0.0874 (downstream) 
  Double_t thick_Ca40 = 0.5065; //
  Double_t thick_Ca48 = 0.565; //
  Double_t thick_Fe54 = 0.04663; // 
  Double_t thick_Ti48 = 0.718;

  //target areal density (g/cm^2)
  Double_t sig_H    = rho_H   * thick_H; 
  Double_t sig_D    = rho_D   * thick_D;
  Double_t sig_Be9  = rho_Be9 * thick_Be9;   // 0.986
  Double_t sig_B10  = rho_B10 * thick_B10;   // 0.576
  Double_t sig_B11  = rho_B11 * thick_B11;   // 0.633
  Double_t sig_C12  = rho_C12 * thick_C12;   // 0.574
  Double_t sig_Al27 = rho_Al27 * thick_Al27; // 0.240 (upstream), 0.236 (downstream)
  Double_t sig_Ca40 = rho_Ca40 * thick_Ca40; // 0.785
  Double_t sig_Ca48 = rho_Ca48 * thick_Ca48; // 1.051
  Double_t sig_Fe54 = rho_Fe54 * thick_Fe54; // 0.367
  Double_t sig_Ti48 = rho_Ti48 * thick_Ti48; // 0.294

  // nuclear transparency factors (prob. that hit proton exits the nucleus)
  Double_t T_H    = 1.;
  Double_t T_D    = 1.; 
  Double_t T_Be9  = 0.6;
  Double_t T_B10  = 0.6;
  Double_t T_B11  = 0.6;
  Double_t T_C12  = 0.6;
  Double_t T_Ca40 = 0.43;
  Double_t T_Ca48 = 0.37;
  Double_t T_Fe54 = 0.36;

  Double_t Transparency(TString target=""){

    // helper function to return nuclear transparency of target
    


    if(target=="h")         return T_H;
    else if(target=="d2")   return T_D;
    else if(target=="Be9")  return T_Be9;
    else if(target=="B10")  return T_B10;
    else if(target=="B11")  return T_B11;
    else if(target=="C12")  return T_C12;
    else if(target=="Ca40") return T_Ca40;
    else if(target=="Ca48") return T_Ca48;
    else if(target=="Fe54") return T_Fe54;
    else return 0.;

  }

  Double_t a2(TString target=""){
    
    // helper function to return A/deuterium ratios for a given target with A nucleons

    // a2 scaling factors (A/d) ratios
    Double_t a2_D    = 1.0; 
    Double_t a2_Be9  = 3.9;
    Double_t a2_B10  = 4.0;
    Double_t a2_B11  = 4.0;
    Double_t a2_C12  = 4.5;
    Double_t a2_Ca40 = 4.5;
    Double_t a2_Ca48 = 4.5;
    Double_t a2_Fe54 = 5.2;

    if(target=="d2")        return a2_D;
    else if(target=="Be9")  return a2_Be9;
    else if(target=="B10")  return a2_B10;
    else if(target=="B11")  return a2_B11;
    else if(target=="C12")  return a2_C12;
    else if(target=="Ca40") return a2_Ca40;
    else if(target=="Ca48") return a2_Ca48;
    else if(target=="Fe54") return a2_Fe54;
    else return 0.;
    
  }

  Double_t sig_A(TString target=""){

    // helper function to return target areal density (g/cm^2)

    if(target=="h")         return sig_H;
    else if(target=="d2")   return sig_D;
    else if(target=="Be9")  return sig_Be9;
    else if(target=="B10")  return sig_B10;
    else if(target=="B11")  return sig_B11;
    else if(target=="C12")  return sig_C12;
    else if(target=="Al27") return sig_Al27;    
    else if(target=="Ca40") return sig_Ca40;
    else if(target=="Ca48") return sig_Ca48;
    else if(target=="Fe54") return sig_Fe54;
    else if(target=="Ti48") return sig_Ti48;
    
    else return 0.;
    
    
  }



  
  // variables to calculate the luminosity data
  Double_t targetfac;
  Double_t luminosity;
  Double_t tgt_areal_density;
  
  // calculate simc cafe production luminosity
  Double_t total_simc_time=0;
  Double_t total_simc_charge=0;
  Double_t total_simc_counts=0;
  Double_t luminosity_simc=0;   

  //SIMC CaFe Rate Estimation for 10.6 GeV beam energy
  // These variables are placeholders for MF, SRC rates from SIMC for each target
  Double_t heep_Ib_simc=0; // beam current [uA], will  be used to scale rates to actual current we get
  Double_t cafe_Ib_simc=0;

  Double_t heep_kin0_rates=0;   // H(e,e'p) elastics kin0 (SHMS angle = 8.3 deg)
  Double_t heep_kin1_rates=0;   // H(e,e'p) elastics kin0 (SHMS angle = 7.5 deg)
  Double_t heep_kin2_rates=0;   // H(e,e'p) elastics kin0 (SHMS angle = 6.8 deg)

  Double_t heep_kin0_time=0;   // H(e,e'p) elastics kin0 (SHMS angle = 8.3 deg)
  Double_t heep_kin1_time=0;   // H(e,e'p) elastics kin0 (SHMS angle = 7.5 deg)
  Double_t heep_kin2_time=0;   // H(e,e'p) elastics kin0 (SHMS angle = 6.8 deg)

  Double_t heep_kin0_counts=0;   // H(e,e'p) elastics kin0 (SHMS angle = 8.3 deg)
  Double_t heep_kin1_counts=0;   // H(e,e'p) elastics kin0 (SHMS angle = 7.5 deg)
  Double_t heep_kin2_counts=0;   // H(e,e'p) elastics kin0 (SHMS angle = 6.8 deg)

  Double_t heep_kin0_charge=0;   // H(e,e'p) elastics kin0 (SHMS angle = 8.3 deg)
  Double_t heep_kin1_charge=0;   // H(e,e'p) elastics kin0 (SHMS angle = 7.5 deg)
  Double_t heep_kin2_charge=0;   // H(e,e'p) elastics kin0 (SHMS angle = 6.8 deg)

  Double_t heep_kin0_lumi_simc=0;   // H(e,e'p) elastics kin0 (SHMS angle = 8.3 deg)
  Double_t heep_kin1_lumi_simc=0;   // H(e,e'p) elastics kin1 (SHMS angle = 7.5 deg)
  Double_t heep_kin2_lumi_simc=0;   // H(e,e'p) elastics kin2 (SHMS angle = 6.8 deg)
  
  
  Double_t simc_cafe_counts=0;    // MF/SRC event rate
  Double_t simc_cafe_rates=0;    // MF/SRC event rate


  //this variable can be for either "prod" or "sample" data replay
  TString replay_type;    // either "prod" (for production) or "sample" for sample replay of 100k events (or any other sample evts)

  
  //Initialization parameters (variables actually used in baseAnalyzer.cpp)
  int run;          // run number
  int evtNum;       // number of events replayed (by the input ROOTfile)
  TString daq_mode;   //"coin" or "singles"
  TString e_arm_name;   // electron arm: "HMS" or "SHMS"
  TString h_arm_name;   //hadron arm
  TString analysis_cut;    // analysis cuts: either "bcm_calib", "lumi", "optics", "heep_singles", "heep_coin", "MF" or "SRC" 
  TString analysis_type;   // analysis prefix for "data" or "simc"  
  Bool_t helicity_flag;     //helicity flag
  TString bcm_type;       // BCM type : "BCM1, BCM2, BCM4A, BCM4B, BCM4C"
  Double_t bcm_thrs;      // BCM current threshold cut (analyze data and scalers ONLY above a certain bcm_thrs, e.g. > 5 uA)

  TString trig_type_single;      // trigger type to actually use when applying pre-scale factor to event weight (only if analyzing singles)
  TString trig_type_coin;      // trigger type to actually use when applying pre-scale factor to event weight   (only if analyzing coincidence)
  
  Bool_t combine_runs_flag;     //flag to combine multiple runs (usually sequential runs @ same kinematics in an experiment) -- not in use currently (but can be added)
  Bool_t skim_flag; 
  TString setting; // this will be an indirect argument (no input by user, but determined by kinematics read from report)
  
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
  TString start_of_run;
  TString end_of_run;

  // Read general info from SIMC input file 
  Double_t tgt_mass_simc;
  Double_t beam_energy_simc;
  Double_t hms_p_simc;
  Double_t hms_angle_simc;
  Double_t shms_p_simc;
  Double_t shms_angle_simc;
  
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
  TString simc_InputFileName_rad;
  TString simc_InputFileName_norad;
  TString simc_ifile;  // simc input file (to read central settings used in simulation)
  
  //Output ROOTfile Name
  TString data_OutputFileName_skim_singles; // only for saving singles skimmed leaf variables (with minimal cuts, like bcm cut and edtm cut) 
  TString data_OutputFileName_skim; // only for saving skimmed leaf variables (with minimal cuts, like bcm cut and edtm cut) 
  TString data_OutputFileName;
  TString simc_OutputFileName_rad;
  TString simc_OutputFileName_norad;
  
  //ROOTfile to store combined hists from different runs
  TString data_OutputFileName_combined; 

  
  //Input parameter controls filenames
  TString main_controls_fname;
  TString input_CutFileName;
  TString input_HBinFileName;
  TString input_FileNamePattern;
  TString input_SIMCinfo_FileName;
  
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

  // counting the events underneath the SHMS Cal E/p (for multi-track eff. correction)
  Double_t single_peak_counts, single_peak_counts_err;
  Double_t multi_peak_counts,  multi_peak_counts_err;
  Double_t multi_track_eff,  multi_track_eff_err;
  
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

  Double_t hdcRes_nbins;
  Double_t hdcRes_xmin;
  Double_t hdcRes_xmax;

  	                  
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

  Double_t pdcRes_nbins;
  Double_t pdcRes_xmin;
  Double_t pdcRes_xmax;	                      

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
  // DATA QUALITY CHECK / CUTS STUDY Histograms
  //----------------------------------------------------------------

  // keep track of total charge
  TH1F *H_total_charge;


  //HMS quality check histos
  TH1F *H_hdcRes_fit[dc_PLANES];
  TH1F *H_hbeta_fit;

  //SHMS quality check histos
  TH1F *H_pdcRes_fit[dc_PLANES];
  TH1F *H_pbeta_fit;
  TH1F *H_pcal_fit;

  // coin time quality check
  TH1F *H_ctime_fit;
  
  // ------ Cuts Quality Check Histos ----
  
  // -- NO CUTS HISTOS --
  // kin
  TH1F *H_ep_ctime_noCUT;
  TH1F *H_the_noCUT;
  TH1F *H_W_noCUT;
  TH1F *H_Q2_noCUT;
  TH1F *H_xbj_noCUT;
  TH1F *H_nu_noCUT;
  TH1F *H_q_noCUT;
  TH1F *H_thq_noCUT;
  TH1F *H_Em_nuc_noCUT;
  TH1F *H_Em_src_noCUT;
  TH1F *H_MM_noCUT;
  TH1F *H_Pm_noCUT;
  TH1F *H_thxq_noCUT;
  TH1F *H_thrq_noCUT;
  TH1F *H_cthrq_noCUT;
  TH1F *H_kf_noCUT;
  TH1F *H_Pf_noCUT;
  TH1F *H_thx_noCUT;

  // recon.
  TH1F *H_eytar_noCUT;  
  TH1F *H_eyptar_noCUT; 
  TH1F *H_exptar_noCUT; 
  TH1F *H_edelta_noCUT;
  TH1F *H_hytar_noCUT;  
  TH1F *H_hyptar_noCUT; 
  TH1F *H_hxptar_noCUT; 
  TH1F *H_hdelta_noCUT;

  // detector
  TH1F *H_pCalEtotTrkNorm_noCUT;
  TH1F *H_pHodBetaTrk_noCUT;
  TH1F *H_pNGCerNpeSum_noCUT;
  TH1F *H_hCalEtotTrkNorm_noCUT;
  TH1F *H_hHodBetaTrk_noCUT;
  TH1F *H_hCerNpeSum_noCUT;  

  // 2d histos
  TH2F *H_hxfp_vs_hyfp_noCUT;
  TH2F *H_exfp_vs_eyfp_noCUT;  
  TH2F *H_hXColl_vs_hYColl_noCUT;
  TH2F *H_eXColl_vs_eYColl_noCUT;
  TH2F *H_Em_nuc_vs_Pm_noCUT;
  TH2F *H_Em_src_vs_Pm_noCUT;
  TH2F *H_Q2_vs_xbj_noCUT;
  TH2F *H_cthrq_vs_Pm_noCUT;
  
  // -- CUTS: ACCEPTANCE CUTS ONLY --
  // kin
  TH1F *H_ep_ctime_ACCP;
  TH1F *H_the_ACCP;
  TH1F *H_W_ACCP;
  TH1F *H_Q2_ACCP;
  TH1F *H_xbj_ACCP;
  TH1F *H_nu_ACCP;
  TH1F *H_q_ACCP;
  TH1F *H_thq_ACCP;
  TH1F *H_Em_nuc_ACCP;
  TH1F *H_Em_src_ACCP;
  TH1F *H_MM_ACCP;
  TH1F *H_Pm_ACCP;
  TH1F *H_thxq_ACCP;
  TH1F *H_thrq_ACCP;
  TH1F *H_cthrq_ACCP; 
  TH1F *H_kf_ACCP;
  TH1F *H_Pf_ACCP;
  TH1F *H_thx_ACCP;

  // recon.
  TH1F *H_eytar_ACCP;  
  TH1F *H_eyptar_ACCP; 
  TH1F *H_exptar_ACCP; 
  TH1F *H_edelta_ACCP;
  TH1F *H_hytar_ACCP;  
  TH1F *H_hyptar_ACCP; 
  TH1F *H_hxptar_ACCP; 
  TH1F *H_hdelta_ACCP;

  // detector
  TH1F *H_pCalEtotTrkNorm_ACCP;
  TH1F *H_pHodBetaTrk_ACCP;
  TH1F *H_pNGCerNpeSum_ACCP;
  TH1F *H_hCalEtotTrkNorm_ACCP;
  TH1F *H_hHodBetaTrk_ACCP;
  TH1F *H_hCerNpeSum_ACCP;  

  // 2d histos
  TH2F *H_hxfp_vs_hyfp_ACCP;
  TH2F *H_exfp_vs_eyfp_ACCP;  
  TH2F *H_hXColl_vs_hYColl_ACCP;
  TH2F *H_eXColl_vs_eYColl_ACCP;
  TH2F *H_Em_nuc_vs_Pm_ACCP;
  TH2F *H_Em_src_vs_Pm_ACCP;
  TH2F *H_Q2_vs_xbj_ACCP;
  TH2F *H_cthrq_vs_Pm_ACCP;
 
    // -- CUTS: ACCEPTANCE + PID CUTS ONLY --
  // kin
  TH1F *H_ep_ctime_ACCP_PID;
  TH1F *H_the_ACCP_PID;
  TH1F *H_W_ACCP_PID;
  TH1F *H_Q2_ACCP_PID;
  TH1F *H_xbj_ACCP_PID;
  TH1F *H_nu_ACCP_PID;
  TH1F *H_q_ACCP_PID;
  TH1F *H_thq_ACCP_PID;
  TH1F *H_Em_nuc_ACCP_PID;
  TH1F *H_Em_src_ACCP_PID;
  TH1F *H_MM_ACCP_PID;
  TH1F *H_Pm_ACCP_PID;
  TH1F *H_thxq_ACCP_PID;
  TH1F *H_thrq_ACCP_PID;
  TH1F *H_cthrq_ACCP_PID; 
  TH1F *H_kf_ACCP_PID;
  TH1F *H_Pf_ACCP_PID;
  TH1F *H_thx_ACCP_PID;

  // recon.
  TH1F *H_eytar_ACCP_PID;  
  TH1F *H_eyptar_ACCP_PID; 
  TH1F *H_exptar_ACCP_PID; 
  TH1F *H_edelta_ACCP_PID;
  TH1F *H_hytar_ACCP_PID;  
  TH1F *H_hyptar_ACCP_PID; 
  TH1F *H_hxptar_ACCP_PID; 
  TH1F *H_hdelta_ACCP_PID;

  // detector
  TH1F *H_pCalEtotTrkNorm_ACCP_PID;
  TH1F *H_pHodBetaTrk_ACCP_PID;
  TH1F *H_pNGCerNpeSum_ACCP_PID;
  TH1F *H_hCalEtotTrkNorm_ACCP_PID;
  TH1F *H_hHodBetaTrk_ACCP_PID;
  TH1F *H_hCerNpeSum_ACCP_PID;  

  // 2d histos
  TH2F *H_hxfp_vs_hyfp_ACCP_PID;
  TH2F *H_exfp_vs_eyfp_ACCP_PID;  
  TH2F *H_hXColl_vs_hYColl_ACCP_PID;
  TH2F *H_eXColl_vs_eYColl_ACCP_PID;
  TH2F *H_Em_nuc_vs_Pm_ACCP_PID;
  TH2F *H_Em_src_vs_Pm_ACCP_PID;
  TH2F *H_Q2_vs_xbj_ACCP_PID;
  TH2F *H_cthrq_vs_Pm_ACCP_PID;
  TH2F *H_ebeta_vs_ctime_ACCP_PID; // C.Y. newly added (feb 14)

  // -- CUTS: ACCEPTANCE + PID CUTS + COIN.TIME ONLY --
  // kin
  TH1F *H_ep_ctime_ACCP_PID_CTIME;
  TH1F *H_the_ACCP_PID_CTIME;
  TH1F *H_W_ACCP_PID_CTIME;
  TH1F *H_Q2_ACCP_PID_CTIME;
  TH1F *H_xbj_ACCP_PID_CTIME;
  TH1F *H_nu_ACCP_PID_CTIME;
  TH1F *H_q_ACCP_PID_CTIME;
  TH1F *H_thq_ACCP_PID_CTIME;
  TH1F *H_Em_nuc_ACCP_PID_CTIME;
  TH1F *H_Em_src_ACCP_PID_CTIME;
  TH1F *H_MM_ACCP_PID_CTIME;
  TH1F *H_Pm_ACCP_PID_CTIME;
  TH1F *H_thxq_ACCP_PID_CTIME;
  TH1F *H_thrq_ACCP_PID_CTIME;
  TH1F *H_cthrq_ACCP_PID_CTIME;  
  TH1F *H_kf_ACCP_PID_CTIME;
  TH1F *H_Pf_ACCP_PID_CTIME;
  TH1F *H_thx_ACCP_PID_CTIME;

  // recon.
  TH1F *H_eytar_ACCP_PID_CTIME;  
  TH1F *H_eyptar_ACCP_PID_CTIME; 
  TH1F *H_exptar_ACCP_PID_CTIME; 
  TH1F *H_edelta_ACCP_PID_CTIME;
  TH1F *H_hytar_ACCP_PID_CTIME;  
  TH1F *H_hyptar_ACCP_PID_CTIME; 
  TH1F *H_hxptar_ACCP_PID_CTIME; 
  TH1F *H_hdelta_ACCP_PID_CTIME;

  // detector
  TH1F *H_pCalEtotTrkNorm_ACCP_PID_CTIME;
  TH1F *H_pHodBetaTrk_ACCP_PID_CTIME;
  TH1F *H_pNGCerNpeSum_ACCP_PID_CTIME;
  TH1F *H_hCalEtotTrkNorm_ACCP_PID_CTIME;
  TH1F *H_hHodBetaTrk_ACCP_PID_CTIME;
  TH1F *H_hCerNpeSum_ACCP_PID_CTIME;  

  // 2d histos
  TH2F *H_hxfp_vs_hyfp_ACCP_PID_CTIME;
  TH2F *H_exfp_vs_eyfp_ACCP_PID_CTIME;  
  TH2F *H_hXColl_vs_hYColl_ACCP_PID_CTIME;
  TH2F *H_eXColl_vs_eYColl_ACCP_PID_CTIME;
  TH2F *H_Em_vs_Pm_ACCP_PID_CTIME;      // newly added (feb 14)
  TH2F *H_Em_nuc_vs_Pm_ACCP_PID_CTIME;
  TH2F *H_Em_src_vs_Pm_ACCP_PID_CTIME;
  TH2F *H_Q2_vs_xbj_ACCP_PID_CTIME;
  TH2F *H_Pm_vs_thrq_ACCP_PID_CTIME;  // newly added (feb14)
  TH2F *H_cthrq_vs_Pm_ACCP_PID_CTIME;

  // -- CUTS: ACCEPTANCE + PID CUTS + COIN.TIME + Q2 CUT ONLY --
  TH1F *H_ep_ctime_ACCP_PID_CTIME_Q2;
  TH1F *H_Q2_ACCP_PID_CTIME_Q2;
  TH1F *H_xbj_ACCP_PID_CTIME_Q2;
  TH1F *H_Em_nuc_ACCP_PID_CTIME_Q2;
  TH1F *H_Em_src_ACCP_PID_CTIME_Q2;
  TH1F *H_Pm_ACCP_PID_CTIME_Q2;
  TH1F *H_thrq_ACCP_PID_CTIME_Q2;
  TH1F *H_cthrq_ACCP_PID_CTIME_Q2;
  TH2F *H_hxfp_vs_hyfp_ACCP_PID_CTIME_Q2;
  TH2F *H_exfp_vs_eyfp_ACCP_PID_CTIME_Q2;  
  TH2F *H_hXColl_vs_hYColl_ACCP_PID_CTIME_Q2;
  TH2F *H_eXColl_vs_eYColl_ACCP_PID_CTIME_Q2;
  TH2F *H_Em_nuc_vs_Pm_ACCP_PID_CTIME_Q2;
  TH2F *H_Em_src_vs_Pm_ACCP_PID_CTIME_Q2;
  TH2F *H_Q2_vs_xbj_ACCP_PID_CTIME_Q2;
  TH2F *H_cthrq_vs_Pm_ACCP_PID_CTIME_Q2;

  //  require MF flag
  // -- CUTS: ACCEPTANCE + PID CUTS + COIN.TIME + Q2 + Em CUT ONLY (MF) --
  TH1F *H_ep_ctime_ACCP_PID_CTIME_Q2_Em;
  TH1F *H_Q2_ACCP_PID_CTIME_Q2_Em;
  TH1F *H_xbj_ACCP_PID_CTIME_Q2_Em;
  TH1F *H_Em_nuc_ACCP_PID_CTIME_Q2_Em;
  TH1F *H_Em_src_ACCP_PID_CTIME_Q2_Em;
  TH1F *H_Pm_ACCP_PID_CTIME_Q2_Em;
  TH1F *H_thrq_ACCP_PID_CTIME_Q2_Em;
  TH1F *H_cthrq_ACCP_PID_CTIME_Q2_Em;
  TH2F *H_hxfp_vs_hyfp_ACCP_PID_CTIME_Q2_Em;
  TH2F *H_exfp_vs_eyfp_ACCP_PID_CTIME_Q2_Em;  
  TH2F *H_hXColl_vs_hYColl_ACCP_PID_CTIME_Q2_Em;
  TH2F *H_eXColl_vs_eYColl_ACCP_PID_CTIME_Q2_Em;
  TH2F *H_Em_nuc_vs_Pm_ACCP_PID_CTIME_Q2_Em;
  TH2F *H_Em_src_vs_Pm_ACCP_PID_CTIME_Q2_Em;
  TH2F *H_Q2_vs_xbj_ACCP_PID_CTIME_Q2_Em;
  TH2F *H_cthrq_vs_Pm_ACCP_PID_CTIME_Q2_Em;
 
  //  require MF flag
  // -- CUTS: ACCEPTANCE + PID CUTS + COIN.TIME + Q2 + Em + Pm CUT ONLY (MF) --
  TH1F *H_ep_ctime_ACCP_PID_CTIME_Q2_Em_Pm;
  TH1F *H_Q2_ACCP_PID_CTIME_Q2_Em_Pm;
  TH1F *H_xbj_ACCP_PID_CTIME_Q2_Em_Pm;
  TH1F *H_Em_nuc_ACCP_PID_CTIME_Q2_Em_Pm;
  TH1F *H_Em_src_ACCP_PID_CTIME_Q2_Em_Pm;
  TH1F *H_Pm_ACCP_PID_CTIME_Q2_Em_Pm;
  TH1F *H_thrq_ACCP_PID_CTIME_Q2_Em_Pm;
  TH1F *H_cthrq_ACCP_PID_CTIME_Q2_Em_Pm; 
  TH2F *H_hxfp_vs_hyfp_ACCP_PID_CTIME_Q2_Em_Pm;
  TH2F *H_exfp_vs_eyfp_ACCP_PID_CTIME_Q2_Em_Pm;  
  TH2F *H_hXColl_vs_hYColl_ACCP_PID_CTIME_Q2_Em_Pm;
  TH2F *H_eXColl_vs_eYColl_ACCP_PID_CTIME_Q2_Em_Pm;
  TH2F *H_Em_nuc_vs_Pm_ACCP_PID_CTIME_Q2_Em_Pm;
  TH2F *H_Em_src_vs_Pm_ACCP_PID_CTIME_Q2_Em_Pm;
  TH2F *H_Q2_vs_xbj_ACCP_PID_CTIME_Q2_Em_Pm;
  TH2F *H_cthrq_vs_Pm_ACCP_PID_CTIME_Q2_Em_Pm;
  
  //  require SRC flag
  // -- CUTS: ACCEPTANCE + PID CUTS + COIN.TIME + Q2 + Xbj CUT ONLY (SRC) --
  TH1F *H_ep_ctime_ACCP_PID_CTIME_Q2_Xbj;
  TH1F *H_Q2_ACCP_PID_CTIME_Q2_Xbj;
  TH1F *H_xbj_ACCP_PID_CTIME_Q2_Xbj;
  TH1F *H_Em_nuc_ACCP_PID_CTIME_Q2_Xbj;
  TH1F *H_Em_src_ACCP_PID_CTIME_Q2_Xbj;
  TH1F *H_Pm_ACCP_PID_CTIME_Q2_Xbj;
  TH1F *H_thrq_ACCP_PID_CTIME_Q2_Xbj;
  TH1F *H_cthrq_ACCP_PID_CTIME_Q2_Xbj;
  TH2F *H_hxfp_vs_hyfp_ACCP_PID_CTIME_Q2_Xbj;
  TH2F *H_exfp_vs_eyfp_ACCP_PID_CTIME_Q2_Xbj;  
  TH2F *H_hXColl_vs_hYColl_ACCP_PID_CTIME_Q2_Xbj;
  TH2F *H_eXColl_vs_eYColl_ACCP_PID_CTIME_Q2_Xbj;
  TH2F *H_Em_nuc_vs_Pm_ACCP_PID_CTIME_Q2_Xbj;
  TH2F *H_Em_src_vs_Pm_ACCP_PID_CTIME_Q2_Xbj;
  TH2F *H_Q2_vs_xbj_ACCP_PID_CTIME_Q2_Xbj;
  TH2F *H_cthrq_vs_Pm_ACCP_PID_CTIME_Q2_Xbj;
 
  //  require SRC flag
  // -- CUTS: ACCEPTANCE + PID CUTS + COIN.TIME + Q2 + Xbj + th_rq CUT ONLY (SRC) --
  TH1F *H_ep_ctime_ACCP_PID_CTIME_Q2_Xbj_thrq;
  TH1F *H_Q2_ACCP_PID_CTIME_Q2_Xbj_thrq;
  TH1F *H_xbj_ACCP_PID_CTIME_Q2_Xbj_thrq;
  TH1F *H_Em_nuc_ACCP_PID_CTIME_Q2_Xbj_thrq;
  TH1F *H_Em_src_ACCP_PID_CTIME_Q2_Xbj_thrq;
  TH1F *H_Pm_ACCP_PID_CTIME_Q2_Xbj_thrq;
  TH1F *H_thrq_ACCP_PID_CTIME_Q2_Xbj_thrq;
  TH1F *H_cthrq_ACCP_PID_CTIME_Q2_Xbj_thrq;
  TH2F *H_hxfp_vs_hyfp_ACCP_PID_CTIME_Q2_Xbj_thrq;
  TH2F *H_exfp_vs_eyfp_ACCP_PID_CTIME_Q2_Xbj_thrq;  
  TH2F *H_hXColl_vs_hYColl_ACCP_PID_CTIME_Q2_Xbj_thrq;
  TH2F *H_eXColl_vs_eYColl_ACCP_PID_CTIME_Q2_Xbj_thrq;
  TH2F *H_Em_nuc_vs_Pm_ACCP_PID_CTIME_Q2_Xbj_thrq;
  TH2F *H_Em_src_vs_Pm_ACCP_PID_CTIME_Q2_Xbj_thrq;
  TH2F *H_Q2_vs_xbj_ACCP_PID_CTIME_Q2_Xbj_thrq;
  TH2F *H_cthrq_vs_Pm_ACCP_PID_CTIME_Q2_Xbj_thrq;
 
   //  require SRC flag
  // -- CUTS: ACCEPTANCE + PID CUTS + COIN.TIME + Q2 + Xbj + th_rq + Pm CUT ONLY (SRC) --
  TH1F *H_ep_ctime_ACCP_PID_CTIME_Q2_Xbj_thrq_Pm;
  TH1F *H_Q2_ACCP_PID_CTIME_Q2_Xbj_thrq_Pm;
  TH1F *H_xbj_ACCP_PID_CTIME_Q2_Xbj_thrq_Pm;
  TH1F *H_Em_nuc_ACCP_PID_CTIME_Q2_Xbj_thrq_Pm;
  TH1F *H_Em_src_ACCP_PID_CTIME_Q2_Xbj_thrq_Pm;
  TH1F *H_Pm_ACCP_PID_CTIME_Q2_Xbj_thrq_Pm;
  TH1F *H_thrq_ACCP_PID_CTIME_Q2_Xbj_thrq_Pm;
  TH1F *H_cthrq_ACCP_PID_CTIME_Q2_Xbj_thrq_Pm; 
  TH2F *H_hxfp_vs_hyfp_ACCP_PID_CTIME_Q2_Xbj_thrq_Pm;
  TH2F *H_exfp_vs_eyfp_ACCP_PID_CTIME_Q2_Xbj_thrq_Pm;  
  TH2F *H_hXColl_vs_hYColl_ACCP_PID_CTIME_Q2_Xbj_thrq_Pm;
  TH2F *H_eXColl_vs_eYColl_ACCP_PID_CTIME_Q2_Xbj_thrq_Pm;
  TH2F *H_Em_nuc_vs_Pm_ACCP_PID_CTIME_Q2_Xbj_thrq_Pm;
  TH2F *H_Em_src_vs_Pm_ACCP_PID_CTIME_Q2_Xbj_thrq_Pm;
  TH2F *H_Q2_vs_xbj_ACCP_PID_CTIME_Q2_Xbj_thrq_Pm;
  TH2F *H_cthrq_vs_Pm_ACCP_PID_CTIME_Q2_Xbj_thrq_Pm;

  //------------------------------------

  //----------------------------------------------------------------
  // Detector Histograms (DATA ONLY): PID / TRACKING EFFICIENCY 
  //----------------------------------------------------------------

  //Coin. Time
  TH1F *H_ep_ctime_total;
  TH1F *H_ep_ctime_real;

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

  //--------------------
  // SHMS specific histos
  // for multiple-track corrections 
  //---------------------
  TH1F *H_multitrack_ep_ctime_notrk_ncut;
  TH1F *H_multitrack_ep_ctime_notrk;      // raw coin time spectrum (no track info)
  // SHMS e/p spectra (with only a cut on the raw coin time w/o track info)
  TH1F *H_multitrack_pCalEtotNorm_full_ncut;     // full E/p spectrum (with only cut on ctime notrk)
  TH1F *H_multitrack_pCalEtotNorm_full;     // full E/p spectrum (with only cut on ctime notrk)
  TH1F *H_multitrack_pCalEtotNorm_peak1;    //  E/p single e- peak selected (with only cut on ctime notrk)
  TH1F *H_multitrack_pCalEtotNorm_multipeaks; //  E/p e- multipeaks(>1) selected (with only cut on ctime notrk)
  
  
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
  TH1F *H_cthrq;
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


  //-- 2D Kinematics Histos --
  TH2F *H_Pm_vs_thrq;
  TH2F *H_Pm_vs_thrq_ps;
  
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

  // 2d missing momentum vs. recoin angle (th_rq)
  TH2F *H_Pm_vs_thrq_rand;
  TH2F *H_Pm_vs_thrq_rand_sub;
  
  //-----------END CREATE HISTOGRAMS-----------



  //----------DATA-RELATED VARIABLES-----------
  TTree *tree;
  TTree *tree_skim;
  TTree *tree_skim_singles;
  Long64_t nentries;
  
  //Set-Up counters for accepted singles triggers
  Double_t total_trig1_singles_accp = 0;
  Double_t total_trig2_singles_accp = 0;
  
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
  Double_t total_edtm_accp_bcm_cut_single = 0;
  Double_t total_trig_accp_bcm_cut_single; //generic acc. trig
  Double_t total_trig_accp_bcm_cut_coin; //generic acc. trig

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

  
  //SHMS singles track eff.
  Double_t p_did_singles = 0;
  Double_t p_should_singles = 0;
  Double_t pTrkEff_singles;
  Double_t pTrkEff_singles_err;
  
  //Computer Live Time 
  Double_t cpuLT_trig_single;       //generic computer live time
  Double_t cpuLT_trig_coin;       //generic computer live time

  Double_t cpuLT_trig_err_Bi_single;  //generic cpu live time error (using binomial statistics)
  Double_t cpuLT_trig_err_Bi_coin;  //generic cpu live time error (using binomial statistics)

  Double_t cpuLT_trig_err_Bay_single;  //generic cpu live time error (using bayesian statistics)
  Double_t cpuLT_trig_err_Bay_coin;  //generic cpu live time error (using bayesian statistics)


  
  Double_t cpuLT_trig1, cpuLT_trig1_err_Bi, cpuLT_trig1_err_Bay;
  Double_t cpuLT_trig2, cpuLT_trig2_err_Bi, cpuLT_trig2_err_Bay;
  Double_t cpuLT_trig3, cpuLT_trig3_err_Bi, cpuLT_trig3_err_Bay;
  Double_t cpuLT_trig4, cpuLT_trig4_err_Bi, cpuLT_trig4_err_Bay;
  Double_t cpuLT_trig5, cpuLT_trig5_err_Bi, cpuLT_trig5_err_Bay;
  Double_t cpuLT_trig6, cpuLT_trig6_err_Bi, cpuLT_trig6_err_Bay;

  //Total Live Time (EDTM)
  Double_t tLT_corr_factor;

  Double_t tLT_trig_single; //generic total live time 
  Double_t tLT_trig_coin; //generic total live time 

  
  Double_t tLT_trig_err_Bi_single;  //generic total live time error (using binomial statistics)
  Double_t tLT_trig_err_Bi_coin;  //generic total live time error (using binomial statistics)

  Double_t tLT_trig_err_Bay_single;  //generic total live time error (using bayesian statistics)
  Double_t tLT_trig_err_Bay_coin;  //generic total live time error (using bayesian statistics)

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

  //generic cuts to enable/disable triggers 1-6 (HMS singles -> htrig, SHMS singles -> ptrig, coin -> ptrig)
  Bool_t c_trig1;    
  Bool_t c_trig2;   
  Bool_t c_trig3;   
  Bool_t c_trig4;   
  Bool_t c_trig5;   
  Bool_t c_trig6;

  Bool_t c_notrig1;    
  Bool_t c_notrig2;   
  Bool_t c_notrig3;   
  Bool_t c_notrig4;   
  Bool_t c_notrig5;   
  Bool_t c_notrig6;
  
  //Pre-Scale factor for each pre-trigger (used in computer/total live time calculation, to account for pre-scaled triggers)
  Float_t Ps_factor_single = 1;  //generic pre-scale factor for singles
  Float_t Ps_factor_coin = 1;  //generic pre-scale factor for coin
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

  //SHMS Heavy Gas Cherenkov (Pi /K separation)
  Bool_t pntrack_cut_flag;
  Bool_t c_pntrack;
  Double_t pntracks;
  
  
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

  // in-plane recoil (undetected) angle, theta_rq [deg]
  Bool_t   thrq_MF_cut_flag;
  Bool_t   c_MF_thrq;
  Double_t c_MF_thrq_min;
  Double_t c_MF_thrq_max;
  
  
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


  // e12-10-003: deuteron d(e,e'p) Kinematic Cuts -----
  Bool_t   Q2_deep_cut_flag;
  Bool_t   c_deep_Q2;
  Double_t c_deep_Q2_min;
  Double_t c_deep_Q2_max;

  Bool_t   Em_deep_cut_flag;
  Bool_t   c_deep_Em;
  Double_t c_deep_Em_min;
  Double_t c_deep_Em_max;
  
  
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


  //-------Collimator Study-------
  Bool_t hmsCollCut_flag;      //flag to enable/disable collimator cut
  Bool_t shmsCollCut_flag;

  Bool_t hmsColl_Cut;
  Bool_t shmsColl_Cut;

  TCutG *hms_Coll_gCut;   //HMS Collimator Graphical Cut
  TCutG *shms_Coll_gCut;  //SHMS Collimator Graphical Cut

  Bool_t hms_coll_cut_bool=0;  //boolean to be saved to skimmed rootfile,so that users may be able to make cut 
  Bool_t shms_coll_cut_bool=0;

  //HMS Octagonal Collimator Size (Each of the octagonal points is a multiple of 1 or 1/2 of these values)
  Double_t hms_hsize = 4.575;  //cm
  Double_t hms_vsize = 11.646;
 
  //SHMS Octagonal Collimator Size (Each of the octagonal points is a multiple of 1 or 1/2 of these values)
  Double_t shms_hsize = 17.;  //cm
  Double_t shms_vsize = 25.;

  //Scaling factor to scale collimator cuts from original size cut
  Double_t hms_scale=1.;   //Default
  Double_t shms_scale=1.;

  //------------------------------

  
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
  Bool_t c_kin_deep_Cuts; //kinematics cuts ( e12-10-003 d(e,e'p) kinematic Cuts )

  
  //------------------END DATA-RELATED VARIABLES DEFINED CUTS-------------------
  
  //-------------TTREE LEAF VARIABLE NAMES (DATA or SIMC)--------------

  //hadron / electron 4-vector components (px,py,pz,E)
  Double_t Pfx, Pfy, Pfz, Ef_p; // hadron
  Double_t kfx, kfy, kfz, Ef_k; // e-
  
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
  Double_t epCoinTime_notrk;
  Double_t eKCoinTime;
  Double_t ePiCoinTime;

  Double_t epCoinTime_center;
  Double_t epCoinTime_center_notrk;

  //HMS DETECTORS
  Double_t hcer_npesum;
  Double_t hcal_etot;
  Double_t hcal_etotnorm;
  Double_t hcal_etottracknorm;
  Double_t hhod_beta_ntrk;
  Double_t hhod_beta;
  Double_t hhod_GoodScinHit;
  Double_t hdc_ntrack;
  Double_t   hdc_TheRealGolden; 
  Double_t hdc_res[dc_PLANES]; 		      		   
  Double_t hdc_nhit[dc_PLANES];		     
 
  //SHMS DETECTORS
  Double_t phgcer_npesum;
  Double_t pngcer_npesum;
  Double_t pcal_etot;
  Double_t pcal_etotnorm;
  Double_t pcal_etottracknorm;
  Double_t phod_beta_ntrk;
  Double_t phod_beta;
  Double_t phod_GoodScinHit;
  Double_t pdc_ntrack;
  Double_t   pdc_TheRealGolden; 
  Double_t pdc_res[dc_PLANES]; 
  Double_t pdc_nhit[dc_PLANES];
 

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
  Double_t ki;
  Double_t kf;       // final electron arm momentum [GeV/c]
  Double_t Pf;       // final proton momentum [GeV/c]
  
  //Additional Kinematics (User-defined)
  Double_t th_x;                   //hadron arm particle central angle
  Double_t MM2;                   //Missing Mass Squared
  Double_t MM_red;               // reduced missing mass (MM - (MA - MP)), A: target nucleus A mass, MP: proton mass

  Double_t Ex;  // energy of detected particle
  Double_t Er;  // energy of recoil system

  
  
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

  Double_t tar_x;  // common tarx in simc
  
  //Collimator Quantities (in spectrometer corrdinate system)
  Double_t hXColl;  //[cm]
  Double_t hYColl;  //[cm]
  Double_t eXColl;  //[cm]
  Double_t eYColl;  //[cm]
  
  //Additional Acceptance Quantities (User-defined)
  Double_t ztar_diff; //[cm]

  //-----------------------------------------------------------------------
  // LOW-LEVEL DRIFT CHAMBER LEAF VARIABLES FOR TRACKING ALGORITHM STUDIES
  //-----------------------------------------------------------------------
  
  Double_t pdc_stubtest;	   
  Double_t pdc_nhits;		   
  Double_t pdc_tnhits;	   
  
  Double_t pdc_chi2dof;
  Int_t    ndata_pdc_track_chisq;
  Double_t pdc_track_chisq[1000];
  Int_t    ndata_pdc_track_nhits;
  Double_t pdc_track_nhits[1000];	   
  Double_t pdc_InsideDipoleExit;
  
  
  Double_t pdc_Ch1_maxhits;	   
  Double_t pdc_Ch1_spacepoints; 
  Double_t pdc_Ch1_nhit;
  Int_t    ndata_pdc_Ch1_ncombos;	   
  Double_t pdc_Ch1_ncombos[1000]; // this has P.dc.Ch1.ncombos[Ndata.P.dc.Ch1.ncombos]	   
  Int_t    ndata_pdc_Ch1_stub_x;
  Double_t pdc_Ch1_stub_x[1000];
  Int_t    ndata_pdc_Ch1_stub_xp; 
  Double_t pdc_Ch1_stub_xp[1000];
  Int_t    ndata_pdc_Ch1_stub_y; 
  Double_t pdc_Ch1_stub_y[1000];
  Int_t    ndata_pdc_Ch1_stub_yp; 
  Double_t pdc_Ch1_stub_yp[1000];	   
  
  Double_t pdc_Ch2_maxhits;	   
  Double_t pdc_Ch2_spacepoints; 
  Double_t pdc_Ch2_nhit;
  Int_t    ndata_pdc_Ch2_ncombos;	   
  Double_t pdc_Ch2_ncombos[1000];
  Int_t    ndata_pdc_Ch2_stub_x; 	   
  Double_t pdc_Ch2_stub_x[1000];
  Int_t    ndata_pdc_Ch2_stub_xp;
  Double_t pdc_Ch2_stub_xp[1000];
  Int_t    ndata_pdc_Ch2_stub_y;	   
  Double_t pdc_Ch2_stub_y[1000];
  Int_t    ndata_pdc_Ch2_stub_yp;
  Double_t pdc_Ch2_stub_yp[1000];	   
  


  
  //static Double_t pdc_rawTDC[dc_PLANES][1000];

  //static Int_t hndata_rawTDC[dc_PLANES];


 

  //----- SIMC Specific TTree Variable Names -----
  Double_t Normfac;
  Double_t Weight;               //This Weight has the cross section in it

  //Thrown quantities (Used to determine spec. resolution)
  Double_t h_deltai;
  Double_t h_yptari;
  Double_t h_xptari;
  Double_t h_ytari;
  
  Double_t e_deltai;
  Double_t e_yptari;
  Double_t e_xptari;
  Double_t e_ytari;
  
  Double_t corrsing;
  Double_t fry;
  Double_t radphot;
  Double_t sigcc;
  Double_t Jacobian;
  Double_t Genweight;
  Double_t SF_weight;
  Double_t Jacobian_corr;
  Double_t sig;
  Double_t sig_recon;
  Double_t sigcc_recon;
  Double_t coul_corr;
  Double_t Ein;                  //single beam energy value (SIMC Uses this energy. If not corr. for energy loss, it should be same as in input file)
  Double_t SF_weight_recon;

  Double_t prob_abs;  // Probability of absorption of particle in the HMS Collimator
                      //(Must be multiplies by the weight. If particle interation is
                      //NOT simulated, it is set to 1.)
  
  //SIMC x-target corrected (used for  X,Y Collimator position calc.) 
  Double_t htarx_corr;
  Double_t etarx_corr;

  
  
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
  Double_t prev_charge_bcm1 = 0.;
  Double_t prev_charge_bcm2 = 0.;
  Double_t prev_charge_bcm4a = 0.;
  Double_t prev_charge_bcm4b = 0.;
  Double_t prev_charge_bcm4c = 0.;
  Double_t prev_s1x_scaler = 0;
  Double_t prev_s1y_scaler = 0;
  Double_t prev_s2x_scaler = 0;
  Double_t prev_s2y_scaler = 0;
  Double_t prev_trig1_scaler = 0;
  Double_t prev_trig2_scaler = 0;
  Double_t prev_trig3_scaler = 0;
  Double_t prev_trig4_scaler = 0;
  Double_t prev_trig5_scaler = 0;
  Double_t prev_trig6_scaler = 0;
  Double_t prev_edtm_scaler = 0;

  //Define Counter Quantities To Store Accumulated Reads
  Double_t total_time = 0.;
  Double_t total_charge_bcm = 0.; // placeholder for arbitrary BCM type (determined by user)
  Double_t total_charge_bcm1 = 0.;
  Double_t total_charge_bcm2 = 0.;
  Double_t total_charge_bcm4a = 0.;
  Double_t total_charge_bcm4b = 0.;
  Double_t total_charge_bcm4c = 0.;
  Double_t total_s1x_scaler = 0;
  Double_t total_s1y_scaler = 0;
  Double_t total_s2x_scaler = 0;
  Double_t total_s2y_scaler = 0;
  Double_t total_trig1_scaler = 0;
  Double_t total_trig2_scaler = 0;
  Double_t total_trig3_scaler = 0;
  Double_t total_trig4_scaler = 0;
  Double_t total_trig5_scaler = 0;
  Double_t total_trig6_scaler = 0;
  Double_t total_edtm_scaler = 0;

  //Store Accumulated Reads if they passed BCM Current Cut
  Double_t total_time_bcm_cut = 0.;
  Double_t total_charge_bcm_cut = 0.;  // placeholder for arbitrary BCM type (determined by user)
  Double_t total_charge_bcm1_cut = 0.;
  Double_t total_charge_bcm2_cut = 0.;
  Double_t total_charge_bcm4a_cut = 0.;
  Double_t total_charge_bcm4b_cut = 0.;
  Double_t total_charge_bcm4c_cut = 0.;
  Double_t total_s1x_scaler_bcm_cut = 0;
  Double_t total_s1y_scaler_bcm_cut = 0;
  Double_t total_s2x_scaler_bcm_cut = 0;
  Double_t total_s2y_scaler_bcm_cut = 0;
  Double_t total_trig1_scaler_bcm_cut = 0;
  Double_t total_trig2_scaler_bcm_cut = 0;
  Double_t total_trig3_scaler_bcm_cut = 0;
  Double_t total_trig4_scaler_bcm_cut = 0;
  Double_t total_trig5_scaler_bcm_cut = 0;
  Double_t total_trig6_scaler_bcm_cut = 0;
  Double_t total_edtm_scaler_bcm_cut = 0;
  Double_t total_trig_scaler_bcm_cut_single; //generic trig scaler count
  Double_t total_trig_scaler_bcm_cut_coin; //generic trig scaler count

  //Store Scaler Rates if current cut passed
  Double_t S1XscalerRate_bcm_cut;
  Double_t S1YscalerRate_bcm_cut;
  Double_t S2XscalerRate_bcm_cut;
  Double_t S2YscalerRate_bcm_cut;

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
  
  Double_t trig_rate_single; //generic single trigger rate 
  Double_t trig_rate_coin; //generic coin trigger rate
  
  //Store Average BCM Current
  Double_t  avg_current_bcm_cut;
  
  //---------END SCALER-RELATED VARIABLES----------


  //--------SCALER TTREE VARIABLE NAMES (DATA)---------

  Double_t Scal_evNum;
  Double_t Scal_BCM_charge;  // generic placeholder for BCM charge (depend on user input)
  Double_t Scal_BCM_current; // generic placeholder for BCM current
  
  // C.Y. Oct 3 : added additional bcm info (to write to report as well)
  Double_t Scal_BCM1_charge; 
  Double_t Scal_BCM1_current;
  Double_t Scal_BCM2_charge; 
  Double_t Scal_BCM2_current;
  Double_t Scal_BCM4A_charge; 
  Double_t Scal_BCM4A_current;
  Double_t Scal_BCM4B_charge; 
  Double_t Scal_BCM4B_current;
  Double_t Scal_BCM4C_charge; 
  Double_t Scal_BCM4C_current;
  
  Double_t Scal_time;
  Double_t S1X_scaler;
  Double_t S1Y_scaler;
  Double_t S2X_scaler;
  Double_t S2Y_scaler;

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
  Double_t PhaseSpace = 1;  //default
  
  
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

  // Quality Check Parameter Values (to be used for storing fit results/write to.csv file)
  
  
  // define peak values (max bin content x-value)
  Double_t ctime_offset_peak_val = 0.0;
  Double_t ctime_offset_peak_notrk_val = 0.0;
  Double_t hms_beta_peak_val = 0.0;
  Double_t shms_beta_peak_val = 0.0;
  Double_t shms_ecal_peak_val = 0.0;
    
  // coin time [ns]
  Double_t ctime_offset;                                                                                                                          
  Double_t ctime_offset_err;
  Double_t ctime_sigma;   
  Double_t ctime_sigma_err;

  // HMS Hodo Beta (track)
  Double_t hbeta_mean;                                                                                                                          
  Double_t hbeta_mean_err;
  Double_t hbeta_sigma;   
  Double_t hbeta_sigma_err;

  // SHMS Hodo Beta (track)
  Double_t pbeta_mean;                                                                                                                          
  Double_t pbeta_mean_err;
  Double_t pbeta_sigma;   
  Double_t pbeta_sigma_err;

  // SHMS cal E_dep / p_track
  Double_t pcal_mean;                                                                                                                          
  Double_t pcal_mean_err;
  Double_t pcal_sigma;   
  Double_t pcal_sigma_err;

  // HMS DC Residuals
  Double_t hdc_res_mean[dc_PLANES];                                                                                                                          
  Double_t hdc_res_mean_err[dc_PLANES];
  Double_t hdc_res_sigma[dc_PLANES];   
  Double_t hdc_res_sigma_err[dc_PLANES];

   // SHMS DC Residuals
  Double_t pdc_res_mean[dc_PLANES];                                                                                                                          
  Double_t pdc_res_mean_err[dc_PLANES];
  Double_t pdc_res_sigma[dc_PLANES];   
  Double_t pdc_res_sigma_err[dc_PLANES];
  
  
  //------VARIABLES USED TO WRITE HISTOGRAMS TO ROOT FILE-------

  //Create Categorical TLists to store histograms based on caterogy
  TList *pid_HList;    //store detector histograms (i.e., coin_time, H_hcerNpeSum, H_eCalEtotNorm, . . .)
  TList *kin_HList;    //store kinematical histograms (i.e., Q2, W, th_e, . . .)
  TList *accp_HList;   //store spectrometer accpetance histograms (focal plane, reconstructed, ztar_diff)

  TList * rand_HList;  // store random coin. background of selected histograms
  TList * randSub_HList;  // store random-subtracted variables of selected histograms

  TList * quality_HList; // store quality-check histos (will NOT be weighted or summed over all runs)
  TList * charge_HList;
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
