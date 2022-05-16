//Script to make comparison between SIMC and Commissioning Data from HallC Spring 2018
//Compare Target Reconstruction/FOCAL PLANE/ Kinematics Variables

void compare_deep(int set, int pm, string model, string rad)
{

  gROOT->SetBatch(kTRUE);  
  gStyle->SetOptStat(1001111);
  
  //Pre-defined SIMC/data root file names containing histogram object to comapare
  //TString simc_filename_fsi = Form("../deep_simc_histos_pm%d_lagetfsi_rad_set%d.root", pm, set);
  //TString simc_filename_pwia = Form("../deep_simc_histos_pm%d_lagetpwia_rad_set%d.root", pm, set);
  TString simc_filename_fsi = "../root_files/pm80_Xsec/deep_simc_histos_pm80_lagetfsi_rad_set1.root";
  TString simc_filename_pwia = "../root_files/pm80_Xsec/deep_simc_histos_pm80_lagetpwia_rad_set1.root ";

  //Data File
  //TString data_filename = Form("../deep_data_histos_pm%d_set%d_combined.root", pm, set); 
  TString data_filename = "../root_files/pm80_Xsec/deep_data_histos_pm80_set1_combined.root";
  
  //Open SIMC/data ROOT files;
  
  TFile *simc_file_fsi = new TFile(simc_filename_fsi);
  TFile *simc_file_pwia = new TFile(simc_filename_pwia);

  TFile *data_file = new TFile(data_filename);

  //---------------Target ----------------
  //Define SIMC histos ('h'-->hadron arm,  'e'-->electron arm)


  //FSI
  TH1F *simc_xtar_fsi =  0;
  TH1F *simc_ytarH_fsi =  0;
  TH1F *simc_ztarH_fsi =  0;
  TH1F *simc_ytarP_fsi =  0;                                                                                                                                     
  TH1F *simc_ztarP_fsi =  0; 

  //PWIA
  TH1F *simc_xtar_pwia =  0;
  TH1F *simc_ytarH_pwia =  0;
  TH1F *simc_ztarH_pwia =  0;
  TH1F *simc_ytarP_pwia =  0;                                                                                                                                     
  TH1F *simc_ztarP_pwia =  0; 



  
  //Define data histos
  TH1F *data_xtarH = 0;
  TH1F *data_ytarH = 0;
  TH1F *data_ztarH = 0;
  
  TH1F *data_xtarP = 0;
  TH1F *data_ytarP = 0;
  TH1F *data_ztarP = 0;
  


  //---------------Target Reconstruction Variables----------------
  //Define SIMC histos ('h'-->hadron arm,  'e'-->electron arm)
  //FSI
  TH1F *simc_eytar_fsi =  0;
  TH1F *simc_exptar_fsi =  0;
  TH1F *simc_eyptar_fsi =  0;
  TH1F *simc_edelta_fsi =  0;

  TH1F *simc_hytar_fsi =  0;
  TH1F *simc_hxptar_fsi =  0;
  TH1F *simc_hyptar_fsi =  0;
  TH1F *simc_hdelta_fsi =  0;

  //PWIA
  TH1F *simc_eytar_pwia =  0;
  TH1F *simc_exptar_pwia =  0;
  TH1F *simc_eyptar_pwia =  0;
  TH1F *simc_edelta_pwia =  0;

  TH1F *simc_hytar_pwia =  0;
  TH1F *simc_hxptar_pwia =  0;
  TH1F *simc_hyptar_pwia =  0;
  TH1F *simc_hdelta_pwia =  0;
  
  //Define data histos
  TH1F *data_eytar = 0;
  TH1F *data_exptar =  0;
  TH1F *data_eyptar =  0;
  TH1F *data_edelta =  0;

  TH1F *data_hytar = 0;
  TH1F *data_hxptar =  0;
  TH1F *data_hyptar =  0;
  TH1F *data_hdelta =  0;

  //-----------------------------------------------------------
 
  //--------------FOCAL PLANE VARIABLES------------------------

 //Define SIMC histos ('h'-->hadron arm,  'e'-->electron arm)
  //FSI
  TH1F *simc_exfp_fsi =  0;
  TH1F *simc_eyfp_fsi =  0;
  TH1F *simc_expfp_fsi =  0;
  TH1F *simc_eypfp_fsi =  0;

  TH1F *simc_hxfp_fsi =  0;
  TH1F *simc_hyfp_fsi =  0;
  TH1F *simc_hxpfp_fsi =  0;
  TH1F *simc_hypfp_fsi =  0;

  //PWIA
  TH1F *simc_exfp_pwia =  0;
  TH1F *simc_eyfp_pwia =  0;
  TH1F *simc_expfp_pwia =  0;
  TH1F *simc_eypfp_pwia =  0;

  TH1F *simc_hxfp_pwia =  0;
  TH1F *simc_hyfp_pwia =  0;
  TH1F *simc_hxpfp_pwia =  0;
  TH1F *simc_hypfp_pwia =  0;
  
  //Define data histos
  TH1F *data_exfp =  0;
  TH1F *data_eyfp =  0;
  TH1F *data_expfp =  0;
  TH1F *data_eypfp =  0;

  TH1F *data_hxfp =  0;
  TH1F *data_hyfp =  0;
  TH1F *data_hxpfp =  0;
  TH1F *data_hypfp =  0;

  //--------------------------------------------------------------

  //-------------------------KINEMATICS---------------------------

  //FSI
  TH1F *simc_Q2_fsi =  0;
  TH1F *simc_omega_fsi =  0;
  TH1F *simc_W2_fsi =  0;
  TH1F *simc_thq_fsi = 0;

  TH1F *simc_xbj_fsi = 0;
  TH1F *simc_th_elec_fsi = 0;                                  
  TH1F *simc_kf_fsi = 0;  
  TH1F *simc_emiss_fsi = 0;

  
  TH1F *simc_Pm_fsi = 0;
  TH1F *simc_Pf_fsi = 0;
  TH1F *simc_th_prot_fsi = 0;
  TH1F *simc_q_fsi = 0;    //q-vector magnitude
  TH1F *simc_thpq_fsi = 0;
  TH1F *simc_thnq_fsi = 0;
  
  
  TH1F *simc_MM_fsi = 0;
  TH1F *simc_En_fsi = 0;
  TH1F *simc_Ep_fsi = 0;
  TH1F *simc_Kn_fsi = 0;
  TH1F *simc_Kp_fsi = 0;
  TH1F *simc_Pmx_fsi = 0;
  TH1F *simc_Pmy_fsi = 0;
  TH1F *simc_Pmz_fsi = 0;

  //PWIA
  TH1F *simc_Q2_pwia =  0;
  TH1F *simc_omega_pwia =  0;
  TH1F *simc_W2_pwia =  0;
  TH1F *simc_thq_pwia = 0;

  TH1F *simc_xbj_pwia = 0;
  TH1F *simc_th_elec_pwia = 0;                                  
  TH1F *simc_kf_pwia = 0;  
  TH1F *simc_emiss_pwia = 0;

  
  TH1F *simc_Pm_pwia = 0;
  TH1F *simc_Pf_pwia = 0;
  TH1F *simc_th_prot_pwia = 0;
  TH1F *simc_q_pwia = 0;    //q-vector magnitude
  TH1F *simc_thpq_pwia = 0;
  TH1F *simc_thnq_pwia = 0;
  
  
  TH1F *simc_MM_pwia = 0;
  TH1F *simc_En_pwia = 0;
  TH1F *simc_Ep_pwia = 0;
  TH1F *simc_Kn_pwia = 0;
  TH1F *simc_Kp_pwia = 0;
  TH1F *simc_Pmx_pwia = 0;
  TH1F *simc_Pmy_pwia = 0;
  TH1F *simc_Pmz_pwia = 0;
  
  
  //Define data histos
  TH1F *data_Q2 =  0;
  TH1F *data_omega =  0;
  TH1F *data_W2 =  0;
  TH1F *data_thq = 0;

  TH1F *data_xbj = 0;
  TH1F *data_th_elec = 0;
  TH1F *data_kf = 0;
  TH1F *data_emiss = 0;

   //Kinematics 2
  TH1F *data_Pm = 0;
  TH1F *data_Pf = 0;
  TH1F *data_th_prot = 0;
  TH1F *data_q = 0;    //q-vector magnitude
  TH1F *data_thpq = 0;
  TH1F *data_thnq = 0;

  //Kinematics 3
  TH1F *data_MM = 0;
  TH1F *data_En = 0;
  TH1F *data_Ep = 0;
  TH1F *data_Kn = 0;
  TH1F *data_Kp = 0;
  TH1F *data_Pmx = 0;
  TH1F *data_Pmy = 0;
  TH1F *data_Pmz = 0;

  //DATA/SIMC Ratios
  TH1F *dataPm_fsi = 0;
  TH1F *dataPm_pwia = 0;

  TH1F *data_thnq_fsi = 0;
  TH1F *data_thnq_pwia = 0;

  //------------Miscellaneous PLots----------

  //FSI
  TH1F *simc_ztar_diff_fsi = 0; 
  TH2F *simc_HMS_Coll_fsi = 0;
  TH2F *simc_SHMS_Coll_fsi = 0;

  //PWIA
  TH1F *simc_ztar_diff_pwia = 0; 
  TH2F *simc_HMS_Coll_pwia = 0;
  TH2F *simc_SHMS_Coll_pwia = 0;
  

  TH1F *data_ztar_diff = 0;
  TH1F *data_CoinTime = 0;
  TH1F *data_pid_eCal = 0;
  TH2F *data_HMS_Coll = 0;
  TH2F *data_SHMS_Coll = 0; 


  
  //---------------------------------------------------------------

  //---------Get Miscellaneous Histograms---------

 //change to FSI simc_file
  simc_file_fsi->cd();

  //Get Histogram objects from SIMC rootfile  
  simc_file_fsi->GetObject("H_ztar_diff", simc_ztar_diff_fsi);
  simc_file_fsi->GetObject("H_hXColl_vs_hYColl", simc_HMS_Coll_fsi);  
  simc_file_fsi->GetObject("H_eXColl_vs_eYColl", simc_SHMS_Coll_fsi); 

  
  //Set SIMC Histo Aesthetics          
  simc_ztar_diff_fsi->SetFillColorAlpha(kRed, 0.35);                                                                                                                                   
  simc_ztar_diff_fsi->SetFillStyle(3004);                                                                                                                                            
  simc_ztar_diff_fsi->SetLineColor(kRed);                                                                                                                                              
                                                                                                                                            
  //change to PWIA simc_file
  simc_file_pwia->cd();

  //Get Histogram objects from SIMC rootfile  
  simc_file_pwia->GetObject("H_ztar_diff", simc_ztar_diff_pwia);
  simc_file_pwia->GetObject("H_hXColl_vs_hYColl", simc_HMS_Coll_pwia);  
  simc_file_pwia->GetObject("H_eXColl_vs_eYColl", simc_SHMS_Coll_pwia); 

  
  //Set SIMC Histo Aesthetics          
  simc_ztar_diff_pwia->SetFillColorAlpha(kBlue, 0.35);                                                                                                                                   
  simc_ztar_diff_pwia->SetFillStyle(3004);                                                                                                                                            
  simc_ztar_diff_pwia->SetLineColor(kBlue);          

  
  //Get Histogram objects from DATA rootfile    
  data_file->GetObject("H_ztar_diff", data_ztar_diff);
  data_file->GetObject("H_hXColl_vs_hYColl", data_HMS_Coll); 
  data_file->GetObject("H_eXColl_vs_eYColl", data_SHMS_Coll); 

  data_file->GetObject("H_pcal_etotTrkNorm", data_pid_eCal); 
  data_file->GetObject("H_ctime", data_CoinTime); 

  //Set DATA Histo Aesthetics                                                                                                                         
  data_ztar_diff->SetLineColor(kBlack);  
  data_ztar_diff->SetLineWidth(2); 

  data_pid_eCal->SetLineColor(kBlack); 
  data_pid_eCal->SetLineWidth(2); 

  data_CoinTime->SetLineColor(kBlack); 
  data_CoinTime->SetLineWidth(2); 


  //----------Get Target Histograms------------------
  //Get Histogram objects from SIMC rootfile

  //===FSI===
  simc_file_fsi->GetObject("H_hx_tar", simc_xtar_fsi);
  simc_file_fsi->GetObject("H_hy_tar", simc_ytarH_fsi);
  simc_file_fsi->GetObject("H_hz_tar", simc_ztarH_fsi);
  simc_file_fsi->GetObject("H_ey_tar", simc_ytarP_fsi);    
  simc_file_fsi->GetObject("H_ez_tar", simc_ztarP_fsi);   

  //Set SIMC Histo Aesthetics
  simc_xtar_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_xtar_fsi->SetFillStyle(3004);
  simc_xtar_fsi->SetLineColor(kRed);

  simc_ytarH_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_ytarH_fsi->SetFillStyle(3004);
  simc_ytarH_fsi->SetLineColor(kRed);

  simc_ztarH_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_ztarH_fsi->SetFillStyle(3004);
  simc_ztarH_fsi->SetLineColor(kRed);

  simc_ytarP_fsi->SetFillColorAlpha(kRed, 0.35);          
  simc_ytarP_fsi->SetFillStyle(3004);                   
  simc_ytarP_fsi->SetLineColor(kRed);
  
  simc_ztarP_fsi->SetFillColorAlpha(kRed, 0.35);                                          
  simc_ztarP_fsi->SetFillStyle(3004); 
  simc_ztarP_fsi->SetLineColor(kRed);

  //===PWIA===
  simc_file_pwia->GetObject("H_hx_tar", simc_xtar_pwia);
  simc_file_pwia->GetObject("H_hy_tar", simc_ytarH_pwia);
  simc_file_pwia->GetObject("H_hz_tar", simc_ztarH_pwia);
  simc_file_pwia->GetObject("H_ey_tar", simc_ytarP_pwia);    
  simc_file_pwia->GetObject("H_ez_tar", simc_ztarP_pwia);   

  //Set SIMC Histo Aesthetics
  simc_xtar_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_xtar_pwia->SetFillStyle(3004);
  simc_xtar_pwia->SetLineColor(kBlue);

  simc_ytarH_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_ytarH_pwia->SetFillStyle(3004);
  simc_ytarH_pwia->SetLineColor(kBlue);

  simc_ztarH_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_ztarH_pwia->SetFillStyle(3004);
  simc_ztarH_pwia->SetLineColor(kBlue);

  simc_ytarP_pwia->SetFillColorAlpha(kBlue, 0.35);          
  simc_ytarP_pwia->SetFillStyle(3004);                   
  simc_ytarP_pwia->SetLineColor(kBlue);
  
  simc_ztarP_pwia->SetFillColorAlpha(kBlue, 0.35);                                          
  simc_ztarP_pwia->SetFillStyle(3004); 
  simc_ztarP_pwia->SetLineColor(kBlue);

  //change to data_file
  data_file->cd();

  //Get Histogram objects from data rootfile
  data_file->GetObject("H_hx_tar", data_xtarH);
  data_file->GetObject("H_hy_tar", data_ytarH);
  data_file->GetObject("H_hz_tar", data_ztarH);

  data_file->GetObject("H_ex_tar", data_xtarP);
  data_file->GetObject("H_ey_tar", data_ytarP);
  data_file->GetObject("H_ez_tar", data_ztarP);


    //Set data Histo Aesthetics
  data_xtarH->SetLineColor(kBlack);
  data_xtarH->SetLineWidth(2);
  data_ytarH->SetLineColor(kBlack);
  data_ytarH->SetLineWidth(2);
  data_ztarH->SetLineColor(kBlack);
  data_ztarH->SetLineWidth(2);


  data_xtarP->SetLineColor(kBlack);
  data_xtarP->SetLineWidth(2);
  data_ytarP->SetLineColor(kBlack);
  data_ytarP->SetLineWidth(2);
  data_ztarP->SetLineColor(kBlack);
  data_ztarP->SetLineWidth(2);

  //-----------------------------------------------------------------


  //---------------------------------------------------------------

  //----------Get Target Reconstructed Histograms------------------

  //===FSI===
 //change to simc_file
  simc_file_fsi->cd();

  //Get Histogram objects from SIMC rootfile
  simc_file_fsi->GetObject("H_eytar", simc_eytar_fsi);
  simc_file_fsi->GetObject("H_exptar", simc_exptar_fsi);
  simc_file_fsi->GetObject("H_eyptar", simc_eyptar_fsi);
  simc_file_fsi->GetObject("H_edelta", simc_edelta_fsi);

  simc_file_fsi->GetObject("H_hytar", simc_hytar_fsi);
  simc_file_fsi->GetObject("H_hxptar", simc_hxptar_fsi);
  simc_file_fsi->GetObject("H_hyptar", simc_hyptar_fsi);
  simc_file_fsi->GetObject("H_hdelta", simc_hdelta_fsi);

  //Set SIMC Histo Aesthetics
  simc_eytar_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_eytar_fsi->SetFillStyle(3004);
  simc_eytar_fsi->SetLineColor(kRed);
  
  simc_exptar_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_exptar_fsi->SetFillStyle(3004);
  simc_exptar_fsi->SetLineColor(kRed);
    
  simc_eyptar_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_eyptar_fsi->SetFillStyle(3004);
  simc_eyptar_fsi->SetLineColor(kRed);

  simc_edelta_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_edelta_fsi->SetFillStyle(3004);
  simc_edelta_fsi->SetLineColor(kRed);

  simc_hytar_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_hytar_fsi->SetFillStyle(3004);
  simc_hytar_fsi->SetLineColor(kRed);

  simc_hxptar_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_hxptar_fsi->SetFillStyle(3004);
  simc_hxptar_fsi->SetLineColor(kRed);

  simc_hyptar_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_hyptar_fsi->SetFillStyle(3004);
  simc_hyptar_fsi->SetLineColor(kRed);

  simc_hdelta_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_hdelta_fsi->SetFillStyle(3004);
  simc_hdelta_fsi->SetLineColor(kRed);

  //===PWIA===
 //change to simc_file
  simc_file_pwia->cd();

  //Get Histogram objects from SIMC rootfile
  simc_file_pwia->GetObject("H_eytar", simc_eytar_pwia);
  simc_file_pwia->GetObject("H_exptar", simc_exptar_pwia);
  simc_file_pwia->GetObject("H_eyptar", simc_eyptar_pwia);
  simc_file_pwia->GetObject("H_edelta", simc_edelta_pwia);

  simc_file_pwia->GetObject("H_hytar", simc_hytar_pwia);
  simc_file_pwia->GetObject("H_hxptar", simc_hxptar_pwia);
  simc_file_pwia->GetObject("H_hyptar", simc_hyptar_pwia);
  simc_file_pwia->GetObject("H_hdelta", simc_hdelta_pwia);

  //Set SIMC Histo Aesthetics
  simc_eytar_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_eytar_pwia->SetFillStyle(3004);
  simc_eytar_pwia->SetLineColor(kBlue);
  
  simc_exptar_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_exptar_pwia->SetFillStyle(3004);
  simc_exptar_pwia->SetLineColor(kBlue);
    
  simc_eyptar_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_eyptar_pwia->SetFillStyle(3004);
  simc_eyptar_pwia->SetLineColor(kBlue);

  simc_edelta_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_edelta_pwia->SetFillStyle(3004);
  simc_edelta_pwia->SetLineColor(kBlue);

  simc_hytar_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_hytar_pwia->SetFillStyle(3004);
  simc_hytar_pwia->SetLineColor(kBlue);

  simc_hxptar_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_hxptar_pwia->SetFillStyle(3004);
  simc_hxptar_pwia->SetLineColor(kBlue);

  simc_hyptar_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_hyptar_pwia->SetFillStyle(3004);
  simc_hyptar_pwia->SetLineColor(kBlue);

  simc_hdelta_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_hdelta_pwia->SetFillStyle(3004);
  simc_hdelta_pwia->SetLineColor(kBlue);


  //change to data_file
  data_file->cd();

  //Get Histogram objects from data rootfile
  data_file->GetObject("H_eytar", data_eytar);
  data_file->GetObject("H_exptar", data_exptar);
  data_file->GetObject("H_eyptar", data_eyptar);
  data_file->GetObject("H_edelta", data_edelta);
  
  data_file->GetObject("H_hytar", data_hytar);
  data_file->GetObject("H_hxptar", data_hxptar);
  data_file->GetObject("H_hyptar", data_hyptar);
  data_file->GetObject("H_hdelta", data_hdelta);

  //Set data Histo Aesthetics
  data_eytar->SetLineColor(kBlack);
  data_eytar->SetLineWidth(2);
  data_exptar->SetLineColor(kBlack);
  data_exptar->SetLineWidth(2);
  data_eyptar->SetLineColor(kBlack);
  data_eyptar->SetLineWidth(2);
  data_edelta->SetLineColor(kBlack);
  data_edelta->SetLineWidth(2);

  data_hytar->SetLineColor(kBlack);
  data_hytar->SetLineWidth(2);
  data_hxptar->SetLineColor(kBlack);
  data_hxptar->SetLineWidth(2);
  data_hyptar->SetLineColor(kBlack);
  data_hyptar->SetLineWidth(2);
  data_hdelta->SetLineColor(kBlack);
  data_hdelta->SetLineWidth(2);

  //-----------------------------------------------------------------

  

  //---------------Get FOCAL PLANE Histograms------------------------

  //===FSI===
  //change to simc_file
  simc_file_fsi->cd();

  //Get Histogram objects from SIMC rootfile
  simc_file_fsi->GetObject("H_exfp", simc_exfp_fsi);
  simc_file_fsi->GetObject("H_eyfp", simc_eyfp_fsi);
  simc_file_fsi->GetObject("H_expfp", simc_expfp_fsi);
  simc_file_fsi->GetObject("H_eypfp", simc_eypfp_fsi);

  simc_file_fsi->GetObject("H_hxfp", simc_hxfp_fsi);
  simc_file_fsi->GetObject("H_hyfp", simc_hyfp_fsi);
  simc_file_fsi->GetObject("H_hxpfp", simc_hxpfp_fsi);
  simc_file_fsi->GetObject("H_hypfp", simc_hypfp_fsi);

  //Set SIMC Histo Aesthetics
  simc_exfp_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_exfp_fsi->SetFillStyle(3004);
  simc_exfp_fsi->SetLineColor(kRed);
  
  simc_eyfp_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_eyfp_fsi->SetFillStyle(3004);
  simc_eyfp_fsi->SetLineColor(kRed);
  
  simc_expfp_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_expfp_fsi->SetFillStyle(3004);
  simc_expfp_fsi->SetLineColor(kRed);

  simc_eypfp_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_eypfp_fsi->SetFillStyle(3004);
  simc_eypfp_fsi->SetLineColor(kRed);

  simc_hxfp_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_hxfp_fsi->SetFillStyle(3004);
  simc_hxfp_fsi->SetLineColor(kRed);

  simc_hyfp_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_hyfp_fsi->SetFillStyle(3004);
  simc_hyfp_fsi->SetLineColor(kRed);

  simc_hxpfp_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_hxpfp_fsi->SetFillStyle(3004);
  simc_hxpfp_fsi->SetLineColor(kRed);

  simc_hypfp_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_hypfp_fsi->SetFillStyle(3004);
  simc_hypfp_fsi->SetLineColor(kRed);

  //===PWIA===
  //change to simc_file
  simc_file_pwia->cd();

  //Get Histogram objects from SIMC rootfile
  simc_file_pwia->GetObject("H_exfp", simc_exfp_pwia);
  simc_file_pwia->GetObject("H_eyfp", simc_eyfp_pwia);
  simc_file_pwia->GetObject("H_expfp", simc_expfp_pwia);
  simc_file_pwia->GetObject("H_eypfp", simc_eypfp_pwia);

  simc_file_pwia->GetObject("H_hxfp", simc_hxfp_pwia);
  simc_file_pwia->GetObject("H_hyfp", simc_hyfp_pwia);
  simc_file_pwia->GetObject("H_hxpfp", simc_hxpfp_pwia);
  simc_file_pwia->GetObject("H_hypfp", simc_hypfp_pwia);

  //Set SIMC Histo Aesthetics
  simc_exfp_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_exfp_pwia->SetFillStyle(3004);
  simc_exfp_pwia->SetLineColor(kBlue);
  
  simc_eyfp_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_eyfp_pwia->SetFillStyle(3004);
  simc_eyfp_pwia->SetLineColor(kBlue);
  
  simc_expfp_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_expfp_pwia->SetFillStyle(3004);
  simc_expfp_pwia->SetLineColor(kBlue);

  simc_eypfp_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_eypfp_pwia->SetFillStyle(3004);
  simc_eypfp_pwia->SetLineColor(kBlue);

  simc_hxfp_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_hxfp_pwia->SetFillStyle(3004);
  simc_hxfp_pwia->SetLineColor(kBlue);

  simc_hyfp_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_hyfp_pwia->SetFillStyle(3004);
  simc_hyfp_pwia->SetLineColor(kBlue);

  simc_hxpfp_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_hxpfp_pwia->SetFillStyle(3004);
  simc_hxpfp_pwia->SetLineColor(kBlue);

  simc_hypfp_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_hypfp_pwia->SetFillStyle(3004);
  simc_hypfp_pwia->SetLineColor(kBlue);

  //change to data_file
  data_file->cd();

  //Get Histogram objects from data rootfile
  data_file->GetObject("H_exfp", data_exfp);
  data_file->GetObject("H_eyfp", data_eyfp);
  data_file->GetObject("H_expfp", data_expfp);
  data_file->GetObject("H_eypfp", data_eypfp);

  data_file->GetObject("H_hxfp", data_hxfp);
  data_file->GetObject("H_hyfp", data_hyfp);
  data_file->GetObject("H_hxpfp", data_hxpfp);
  data_file->GetObject("H_hypfp", data_hypfp);
  //Set data Histo Aesthetics
  data_exfp->SetLineColor(kBlack);
  data_exfp->SetLineWidth(2);
  data_eyfp->SetLineColor(kBlack);
  data_eyfp->SetLineWidth(2);
  data_expfp->SetLineColor(kBlack);
  data_expfp->SetLineWidth(2);
  data_eypfp->SetLineColor(kBlack);
  data_eypfp->SetLineWidth(2);

  data_hxfp->SetLineColor(kBlack);
  data_hxfp->SetLineWidth(2);
  data_hyfp->SetLineColor(kBlack);
  data_hyfp->SetLineWidth(2);
  data_hxpfp->SetLineColor(kBlack);
  data_hxpfp->SetLineWidth(2);
  data_hypfp->SetLineColor(kBlack);
  data_hypfp->SetLineWidth(2);

  //--------------------------------------------------------------
  
  //------------------Get KINEMATICS VARIABLES--------------------

  //===FSI===
  //change to simc_file
  simc_file_fsi->cd();

  //Get Histogram objects from SIMC rootfile
  simc_file_fsi->GetObject("H_Q2", simc_Q2_fsi);
  simc_file_fsi->GetObject("H_omega", simc_omega_fsi);
  simc_file_fsi->GetObject("H_W2", simc_W2_fsi);
  simc_file_fsi->GetObject("H_theta_q", simc_thq_fsi);

  simc_file_fsi->GetObject("H_xbj", simc_xbj_fsi);
  simc_file_fsi->GetObject("H_theta_elec", simc_th_elec_fsi);
  simc_file_fsi->GetObject("H_kf", simc_kf_fsi);
  simc_file_fsi->GetObject("H_Emiss", simc_emiss_fsi);

  simc_file_fsi->GetObject("H_Pm", simc_Pm_fsi);
  simc_file_fsi->GetObject("H_Pf", simc_Pf_fsi);
  simc_file_fsi->GetObject("H_theta_prot", simc_th_prot_fsi);
  simc_file_fsi->GetObject("H_q", simc_q_fsi);
  simc_file_fsi->GetObject("H_theta_pq", simc_thpq_fsi);
  simc_file_fsi->GetObject("H_theta_nq", simc_thnq_fsi);
  
  simc_file_fsi->GetObject("H_MM", simc_MM_fsi);
  simc_file_fsi->GetObject("H_En", simc_En_fsi);
  simc_file_fsi->GetObject("H_Ep", simc_Ep_fsi);
  simc_file_fsi->GetObject("H_Kn", simc_Kn_fsi);
  simc_file_fsi->GetObject("H_Kp", simc_Kp_fsi);
  simc_file_fsi->GetObject("H_Pmx_Lab", simc_Pmx_fsi);
  simc_file_fsi->GetObject("H_Pmy_Lab", simc_Pmy_fsi);
  simc_file_fsi->GetObject("H_Pmz_Lab", simc_Pmz_fsi);

  //Set SIMC Histo Aesthetics
  simc_Q2_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_Q2_fsi->SetFillStyle(3004);
  simc_Q2_fsi->SetLineColor(kRed);
  
  simc_omega_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_omega_fsi->SetFillStyle(3004);
  simc_omega_fsi->SetLineColor(kRed);
  
  simc_W2_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_W2_fsi->SetFillStyle(3004);
  simc_W2_fsi->SetLineColor(kRed);

  simc_thq_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_thq_fsi->SetFillStyle(3004);
  simc_thq_fsi->SetLineColor(kRed);

  simc_xbj_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_xbj_fsi->SetFillStyle(3004);
  simc_xbj_fsi->SetLineColor(kRed);

  simc_th_elec_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_th_elec_fsi->SetFillStyle(3004);
  simc_th_elec_fsi->SetLineColor(kRed);

  simc_kf_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_kf_fsi->SetFillStyle(3004);
  simc_kf_fsi->SetLineColor(kRed);
    
  simc_emiss_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_emiss_fsi->SetFillStyle(3004);
  simc_emiss_fsi->SetLineColor(kRed);

  simc_Pm_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_Pm_fsi->SetFillStyle(3004);
  simc_Pm_fsi->SetLineColor(kRed);

  simc_Pf_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_Pf_fsi->SetFillStyle(3004);
  simc_Pf_fsi->SetLineColor(kRed);

  simc_th_prot_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_th_prot_fsi->SetFillStyle(3004);
  simc_th_prot_fsi->SetLineColor(kRed);

  simc_q_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_q_fsi->SetFillStyle(3004);
  simc_q_fsi->SetLineColor(kRed);

  simc_thpq_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_thpq_fsi->SetFillStyle(3004);
  simc_thpq_fsi->SetLineColor(kRed);

  simc_thnq_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_thnq_fsi->SetFillStyle(3004);
  simc_thnq_fsi->SetLineColor(kRed);

  simc_MM_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_MM_fsi->SetFillStyle(3004);
  simc_MM_fsi->SetLineColor(kRed);

  simc_En_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_En_fsi->SetFillStyle(3004);
  simc_En_fsi->SetLineColor(kRed);

  simc_Ep_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_Ep_fsi->SetFillStyle(3004);
  simc_Ep_fsi->SetLineColor(kRed);

  simc_Kn_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_Kn_fsi->SetFillStyle(3004);
  simc_Kn_fsi->SetLineColor(kRed);

  simc_Kp_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_Kp_fsi->SetFillStyle(3004);
  simc_Kp_fsi->SetLineColor(kRed);

  simc_Pmx_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_Pmx_fsi->SetFillStyle(3004);
  simc_Pmx_fsi->SetLineColor(kRed);

  simc_Pmy_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_Pmy_fsi->SetFillStyle(3004);
  simc_Pmy_fsi->SetLineColor(kRed);

  simc_Pmz_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_Pmz_fsi->SetFillStyle(3004);
  simc_Pmz_fsi->SetLineColor(kRed);

  //===PWIA===
  //change to simc_file
  simc_file_pwia->cd();

  //Get Histogram objects from SIMC rootfile
  simc_file_pwia->GetObject("H_Q2", simc_Q2_pwia);
  simc_file_pwia->GetObject("H_omega", simc_omega_pwia);
  simc_file_pwia->GetObject("H_W2", simc_W2_pwia);
  simc_file_pwia->GetObject("H_theta_q", simc_thq_pwia);

  simc_file_pwia->GetObject("H_xbj", simc_xbj_pwia);
  simc_file_pwia->GetObject("H_theta_elec", simc_th_elec_pwia);
  simc_file_pwia->GetObject("H_kf", simc_kf_pwia);
  simc_file_pwia->GetObject("H_Emiss", simc_emiss_pwia);

  simc_file_pwia->GetObject("H_Pm", simc_Pm_pwia);
  simc_file_pwia->GetObject("H_Pf", simc_Pf_pwia);
  simc_file_pwia->GetObject("H_theta_prot", simc_th_prot_pwia);
  simc_file_pwia->GetObject("H_q", simc_q_pwia);
  simc_file_pwia->GetObject("H_theta_pq", simc_thpq_pwia);
  simc_file_pwia->GetObject("H_theta_nq", simc_thnq_pwia);
  
  simc_file_pwia->GetObject("H_MM", simc_MM_pwia);
  simc_file_pwia->GetObject("H_En", simc_En_pwia);
  simc_file_pwia->GetObject("H_Ep", simc_Ep_pwia);
  simc_file_pwia->GetObject("H_Kn", simc_Kn_pwia);
  simc_file_pwia->GetObject("H_Kp", simc_Kp_pwia);
  simc_file_pwia->GetObject("H_Pmx_Lab", simc_Pmx_pwia);
  simc_file_pwia->GetObject("H_Pmy_Lab", simc_Pmy_pwia);
  simc_file_pwia->GetObject("H_Pmz_Lab", simc_Pmz_pwia);

  //Set SIMC Histo Aesthetics
  simc_Q2_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_Q2_pwia->SetFillStyle(3004);
  simc_Q2_pwia->SetLineColor(kBlue);
  
  simc_omega_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_omega_pwia->SetFillStyle(3004);
  simc_omega_pwia->SetLineColor(kBlue);
  
  simc_W2_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_W2_pwia->SetFillStyle(3004);
  simc_W2_pwia->SetLineColor(kBlue);

  simc_thq_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_thq_pwia->SetFillStyle(3004);
  simc_thq_pwia->SetLineColor(kBlue);

  simc_xbj_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_xbj_pwia->SetFillStyle(3004);
  simc_xbj_pwia->SetLineColor(kBlue);

  simc_th_elec_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_th_elec_pwia->SetFillStyle(3004);
  simc_th_elec_pwia->SetLineColor(kBlue);

  simc_kf_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_kf_pwia->SetFillStyle(3004);
  simc_kf_pwia->SetLineColor(kBlue);
    
  simc_emiss_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_emiss_pwia->SetFillStyle(3004);
  simc_emiss_pwia->SetLineColor(kBlue);

  simc_Pm_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_Pm_pwia->SetFillStyle(3004);
  simc_Pm_pwia->SetLineColor(kBlue);

  simc_Pf_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_Pf_pwia->SetFillStyle(3004);
  simc_Pf_pwia->SetLineColor(kBlue);

  simc_th_prot_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_th_prot_pwia->SetFillStyle(3004);
  simc_th_prot_pwia->SetLineColor(kBlue);

  simc_q_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_q_pwia->SetFillStyle(3004);
  simc_q_pwia->SetLineColor(kBlue);

  simc_thpq_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_thpq_pwia->SetFillStyle(3004);
  simc_thpq_pwia->SetLineColor(kBlue);

  simc_thnq_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_thnq_pwia->SetFillStyle(3004);
  simc_thnq_pwia->SetLineColor(kBlue);

  simc_MM_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_MM_pwia->SetFillStyle(3004);
  simc_MM_pwia->SetLineColor(kBlue);

  simc_En_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_En_pwia->SetFillStyle(3004);
  simc_En_pwia->SetLineColor(kBlue);

  simc_Ep_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_Ep_pwia->SetFillStyle(3004);
  simc_Ep_pwia->SetLineColor(kBlue);

  simc_Kn_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_Kn_pwia->SetFillStyle(3004);
  simc_Kn_pwia->SetLineColor(kBlue);

  simc_Kp_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_Kp_pwia->SetFillStyle(3004);
  simc_Kp_pwia->SetLineColor(kBlue);

  simc_Pmx_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_Pmx_pwia->SetFillStyle(3004);
  simc_Pmx_pwia->SetLineColor(kBlue);

  simc_Pmy_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_Pmy_pwia->SetFillStyle(3004);
  simc_Pmy_pwia->SetLineColor(kBlue);

  simc_Pmz_pwia->SetFillColorAlpha(kBlue, 0.35);
  simc_Pmz_pwia->SetFillStyle(3004);
  simc_Pmz_pwia->SetLineColor(kBlue);

  //change to data_file
  data_file->cd();
  
  //Get Histogram objects from data rootfile
  data_file->GetObject("H_Q2", data_Q2);
  data_file->GetObject("H_omega", data_omega);
  data_file->GetObject("H_W2", data_W2);
  data_file->GetObject("H_theta_q", data_thq);

  
  data_file->GetObject("H_xbj", data_xbj);
  data_file->GetObject("H_theta_elec", data_th_elec);
  data_file->GetObject("H_kf", data_kf);
  data_file->GetObject("H_Em_nuc", data_emiss);

  data_file->GetObject("H_Pm", data_Pm);
  data_file->GetObject("H_Pf", data_Pf);
  data_file->GetObject("H_theta_prot", data_th_prot);
  data_file->GetObject("H_q", data_q);
  data_file->GetObject("H_theta_pq", data_thpq);
  data_file->GetObject("H_theta_nq", data_thnq);

  data_file->GetObject("H_MM", data_MM);
  data_file->GetObject("H_En", data_En);
  data_file->GetObject("H_Ep", data_Ep);
  data_file->GetObject("H_Kn", data_Kn);
  data_file->GetObject("H_Kp", data_Kp);
  data_file->GetObject("H_Pmx_Lab", data_Pmx);
  data_file->GetObject("H_Pmy_Lab", data_Pmy);
  data_file->GetObject("H_Pmz_Lab", data_Pmz);

  //Set data Histo Aesthetics
  data_Q2->SetLineColor(kBlack);
  data_Q2->SetLineWidth(2);
  data_omega->SetLineColor(kBlack);
  data_omega->SetLineWidth(2);
  data_W2->SetLineColor(kBlack);
  data_W2->SetLineWidth(2);
  data_thq->SetLineColor(kBlack);
  data_thq->SetLineWidth(2);

  data_xbj->SetLineColor(kBlack);
  data_xbj->SetLineWidth(2);
  data_th_elec->SetLineColor(kBlack);
  data_th_elec->SetLineWidth(2);
  data_kf->SetLineColor(kBlack);
  data_kf->SetLineWidth(2);
  data_emiss->SetLineColor(kBlack);
  data_emiss->SetLineWidth(2);

  data_Pm->SetLineColor(kBlack);
  data_Pm->SetLineWidth(2);
  data_Pf->SetLineColor(kBlack);
  data_Pf->SetLineWidth(2);
  data_th_prot->SetLineColor(kBlack);
  data_th_prot->SetLineWidth(2);
  data_q->SetLineColor(kBlack);
  data_q->SetLineWidth(2);
  data_thpq->SetLineColor(kBlack);
  data_thpq->SetLineWidth(2);
  data_thnq->SetLineColor(kBlack);
  data_thnq->SetLineWidth(2);

  data_MM->SetLineColor(kBlack);
  data_MM->SetLineWidth(2);
  data_En->SetLineColor(kBlack);
  data_En->SetLineWidth(2);
  data_Ep->SetLineColor(kBlack);
  data_Ep->SetLineWidth(2);
  data_Kn->SetLineColor(kBlack);
  data_Kn->SetLineWidth(2);
  data_Kp->SetLineColor(kBlack);
  data_Kp->SetLineWidth(2);
  data_Pmx->SetLineColor(kBlack);
  data_Pmx->SetLineWidth(2);
  data_Pmy->SetLineColor(kBlack);
  data_Pmy->SetLineWidth(2);
  data_Pmz->SetLineColor(kBlack);
  data_Pmz->SetLineWidth(2);

  //Ratios
  data_file->GetObject("H_Pm", dataPm_fsi);
  data_file->GetObject("H_Pm", dataPm_pwia);
  data_file->GetObject("H_theta_nq", data_thnq_fsi);
  data_file->GetObject("H_theta_nq", data_thnq_pwia);

  //Overlay SIMC/data plots (*** VERY IMPORTANT ***: Range and #bins must be same)

  //Plot Miscellaneous
  auto leg_hColl_fsi = new TLegend(0.1,0.8,0.28,0.9);                                                                                             auto leg_hColl_pwia = new TLegend(0.1,0.8,0.28,0.9);  
  auto leg_ztdiff = new TLegend(0.1,0.8,0.28,0.9);        
  auto leg_eCal = new TLegend(0.1,0.8,0.28,0.9);                                                                    
  auto leg_ctime = new TLegend(0.1,0.8,0.28,0.9); 

  TCanvas *c_coll = new TCanvas("c_coll", "Collimator", 5000, 3000);                    
                                    
  c_coll->Divide(2,2);                                                
  c_coll->cd(1);                                                            
  data_HMS_Coll->Draw("colz");             
  c_coll->cd(2); 
  data_SHMS_Coll->Draw("colz"); 
  c_coll->cd(3);
  simc_HMS_Coll_fsi->Draw("colz");  
  leg_hColl_fsi->AddEntry(simc_HMS_Coll_fsi, "SIMC: FSI", "f");
  leg_hColl_fsi->Draw();
  c_coll->cd(4);
  simc_SHMS_Coll_pwia->Draw("colz"); 
  leg_hColl_pwia->AddEntry(simc_HMS_Coll_pwia, "SIMC: PWIA", "f");    
  leg_hColl_pwia->Draw();
  c_coll->SaveAs(Form("Collimator_pm%d_set%d.pdf", pm, set)); 

  TCanvas *c_pid = new TCanvas("c_pid", "PID", 5000, 3000);
  c_pid->Divide(1,3);
  
  c_pid->cd(1);
  simc_ztar_diff_fsi->Draw("histE0");
  data_ztar_diff->Draw("samesE0"); 
  leg_ztdiff->AddEntry(data_ztar_diff, "Data", "f");
  leg_ztdiff->AddEntry(simc_ztar_diff_fsi, "SIMC: FSI", "f");
  simc_ztar_diff_pwia->Draw("sameshistE0");
  leg_ztdiff->AddEntry(simc_ztar_diff_pwia, "SIMC: PWIA", "f");

  leg_ztdiff->Draw();
 
  c_pid->cd(2); 
  data_CoinTime->Draw("E0");
  leg_ctime->AddEntry(data_CoinTime, "Data", "f");        
  leg_ctime->Draw(); 
  
  c_pid->cd(3);                                                                                                                
  data_pid_eCal->Draw("E0");                                                                                                         
  leg_eCal->AddEntry(data_pid_eCal, "Data", "f");                                                                                                                                           
  leg_eCal->Draw();

  c_pid->SaveAs(Form("PID_pm%d_set%d.pdf", pm, set));   
                                                                                                                                                                        
                                                                                                                                                                                                                                                                                           

   //Set Legend
   auto leg5 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg6 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg7 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg8 = new TLegend(0.1,0.8,0.28,0.9);

  
   //-----------------PLOT Target Reconstructed Variables SIMC/Data comparison-----------------------

   //Create A Canvas to store Target Recon. variable comparisons in HADRON ARM
   
   TCanvas *c1 = new TCanvas("c1", "Electron Arm: Target Reconstruction", 5000, 3000);
   c1->Divide(2,2);

   c1->cd(1);
   simc_eytar_fsi->Draw("histE0");
   data_eytar->Draw("samesE0");
   leg5->AddEntry(data_eytar,"Data","f");
   leg5->AddEntry(simc_eytar_fsi,"SIMC: FSI");
   simc_eytar_pwia->Draw("sameshistE0");
   leg5->AddEntry(simc_eytar_pwia,"SIMC: PWIA");
   leg5->Draw();

   c1->cd(2);
   simc_exptar_fsi->Draw("histE0");
   data_exptar->Draw("samesE0");
   leg5->AddEntry(data_exptar,"Data", "f");
   leg5->AddEntry(simc_exptar_fsi,"SIMC: FSI");
   simc_exptar_pwia->Draw("sameshistE0");
   leg5->AddEntry(simc_exptar_pwia,"SIMC: PWIA");
   leg5->Draw();

   c1->cd(3);
   simc_eyptar_fsi->Draw("histE0");
   data_eyptar->Draw("samesE0");
   leg7->AddEntry(data_eyptar,"Data", "f");
   leg7->AddEntry(simc_eyptar_fsi,"SIMC: FSI");
   simc_eyptar_pwia->Draw("sameshistE0");
   leg7->AddEntry(simc_eyptar_pwia,"SIMC: PWIA");
   leg7->Draw();
     
   c1->cd(4);
   simc_edelta_fsi->Draw("histE0");
   data_edelta->Draw("samesE0");
   leg8->AddEntry(data_edelta,"Data", "f");
   leg8->AddEntry(simc_edelta_fsi,"SIMC: FSI");
   simc_edelta_pwia->Draw("sameshistE0");
   leg8->AddEntry(simc_edelta_pwia,"SIMC: PWIA");
   leg8->Draw();

   c1->SaveAs(Form("eArm_TargRecon_pm%d_set%d.pdf", pm, set));

   //------------------------------------------------------------------------------

   
   //-----------------PLOT FOCAL PLANE  Variables SIMC/Data comparison-----------------------

  //Set Legend
   auto leg9 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg10 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg11 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg12 = new TLegend(0.1,0.8,0.28,0.9);

      //Set Legend
   auto leg13 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg14 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg15 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg16 = new TLegend(0.1,0.8,0.28,0.9);

   TCanvas *c2 = new TCanvas("c2", "Electron Arm: Focal Plane", 5000, 3000);
   c2->Divide(2,2);

   c2->cd(1);
   simc_exfp_fsi->Draw("histE0");
   data_exfp->Draw("samesE0");
   leg13->AddEntry(data_exfp,"Data","f");
   leg13->AddEntry(simc_exfp_fsi,"SIMC: FSI");
   simc_exfp_pwia->Draw("sameshistE0");
   leg13->AddEntry(simc_exfp_pwia,"SIMC: PWIA");
   leg13->Draw();
   
   c2->cd(2);
   simc_eyfp_fsi->Draw("histE0");
   data_eyfp->Draw("samesE0");
   leg14->AddEntry(data_eyfp,"Data", "f");
   leg14->AddEntry(simc_eyfp_fsi,"SIMC: FSI");
   simc_eyfp_pwia->Draw("sameshistE0");
   leg14->AddEntry(simc_eyfp_pwia,"SIMC: PWIA");
   leg14->Draw();

   c2->cd(3);
   simc_expfp_fsi->Draw("histE0");
   data_expfp->Draw("samesE0");
   leg15->AddEntry(data_expfp,"Data", "f");
   leg15->AddEntry(simc_expfp_fsi,"SIMC: FSI");
   simc_expfp_pwia->Draw("sameshistE0");
   leg15->AddEntry(simc_expfp_pwia,"SIMC: PWIA");
   leg15->Draw();
     
   c2->cd(4);
   simc_eypfp_fsi->Draw("histE0");
   data_eypfp->Draw("samesE0");
   leg16->AddEntry(data_eypfp,"Data", "f");
   leg16->AddEntry(simc_eypfp_fsi,"SIMC: FSI");
   simc_eypfp_pwia->Draw("sameshistE0");
   leg16->AddEntry(simc_eypfp_pwia,"SIMC: PWIA");
   leg16->Draw();

   c2->SaveAs(Form("eArm_FocalPlane_pm%d_set%d.pdf", pm, set));                                                                                   

   //----------------------------------------------------------- 
 
   
   //-----------------PLOT KINEMATICS SIMC/Data comparison---------------

   //Kinematics 1:  Missing Varibales
   auto leg_em = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_MM = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_Pm = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_Pmx = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_Pmy = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_Pmz = new TLegend(0.1,0.8,0.28,0.9);


   TCanvas *ck1 = new TCanvas("ck1", "Kinematics1", 5000, 3000);
   ck1->Divide(3,2);
   ck1->cd(1);
   simc_emiss_fsi->Draw("histE0");
   data_emiss->Draw("samesE0");
   data_emiss->GetXaxis()->SetTitle("Missing Energy, E_{m} [GeV] ");
   data_emiss->GetXaxis()->CenterTitle();      
   leg_em->AddEntry(data_emiss,"Data","f");
   leg_em->AddEntry(simc_emiss_fsi,"SIMC: FSI");
   simc_emiss_pwia->Draw("sameshistE0");
   leg_em->AddEntry(simc_emiss_pwia,"SIMC: PWIA");
   leg_em->Draw();
   
   ck1->cd(2);
   data_MM->GetXaxis()->SetTitle("Missing Mass, M_{miss} [GeV]");
   data_MM->GetXaxis()->CenterTitle();
   simc_MM_fsi->Draw("histE0");
   data_MM->Draw("samesE0");
   leg_MM->AddEntry(data_MM,"Data", "f");
   leg_MM->AddEntry(simc_MM_fsi,"SIMC: FSI");
   simc_MM_pwia->Draw("sameshistE0");
   leg_MM->AddEntry(simc_MM_pwia,"SIMC: PWIA");
   leg_MM->Draw();

   ck1->cd(3);
   data_Pm->GetXaxis()->SetTitle("Missing Momentum, P_{miss} [GeV]");
   data_Pm->GetXaxis()->CenterTitle();
   simc_Pm_fsi->Draw("histE0");
   data_Pm->Draw("samesE0");
   leg_Pm->AddEntry(data_Pm,"Data", "f");
   leg_Pm->AddEntry(simc_Pm_fsi,"SIMC: FSI");
   simc_Pm_pwia->Draw("sameshistE0");
   leg_Pm->AddEntry(simc_Pm_pwia,"SIMC: PWIA");
   leg_Pm->Draw();


   ck1->cd(4);
   data_Pmx->GetXaxis()->SetTitle("Missing Momentum X-comp., Pm_{x} [GeV]");
   data_Pmx->GetXaxis()->CenterTitle();
   simc_Pmx_fsi->Draw("histE0");
   data_Pmx->Draw("samesE0");
   leg_Pmx->AddEntry(data_Pmx,"Data", "f");
   leg_Pmx->AddEntry(simc_Pmx_fsi,"SIMC: FSI");
   simc_Pmx_pwia->Draw("sameshistE0");
   leg_Pmx->AddEntry(simc_Pmx_pwia,"SIMC: PWIA");
   leg_Pmx->Draw();
   
   ck1->cd(5);
   data_Pmy->GetXaxis()->SetTitle("Missing Momentum Y-comp., Pm_{y} [GeV]");
   data_Pmy->GetXaxis()->CenterTitle();
   simc_Pmy_fsi->Draw("histE0");
   data_Pmy->Draw("samesE0");
   leg_Pmy->AddEntry(data_Pmy,"Data", "f");
   leg_Pmy->AddEntry(simc_Pmy_fsi,"SIMC: FSI");
   simc_Pmy_pwia->Draw("sameshistE0");
   leg_Pmy->AddEntry(simc_Pmy_pwia,"SIMC: PWIA");
   leg_Pmy->Draw();

   ck1->cd(6);
   data_Pmz->GetXaxis()->SetTitle("Missing Momentum Z-comp., Pm_{z} [GeV]");
   data_Pmz->GetXaxis()->CenterTitle();
   simc_Pmz_fsi->Draw("histE0");
   data_Pmz->Draw("samesE0");
   leg_Pmz->AddEntry(data_Pmz,"Data", "f");
   leg_Pmz->AddEntry(simc_Pmz_fsi,"SIMC: FSI");
   simc_Pmz_pwia->Draw("sameshistE0");
   leg_Pmz->AddEntry(simc_Pmz_pwia,"SIMC: PWIA");
   leg_Pmz->Draw();

   ck1->SaveAs(Form("Kinematics1_pm%d_set%d.pdf", pm, set));                                                                   


   //Kinematics 2:  Electron Kinematics

   //Set Legend
   auto leg_Q2 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_om = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_xbj = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_W2 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_the = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_kf = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_thq = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_q = new TLegend(0.1,0.8,0.28,0.9);


   TCanvas *ck2 = new TCanvas("ck2", "Kinematics2", 5000, 3000);
   ck2->Divide(4,2);
   
   ck2->cd(1);
   data_Q2->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
   data_Q2->GetXaxis()->CenterTitle();
   simc_Q2_fsi->Draw("histE0");
   data_Q2->Draw("samesE0");
   leg_Q2->AddEntry(data_Q2,"Data", "f");
   leg_Q2->AddEntry(simc_Q2_fsi,"SIMC: FSI");
   simc_Q2_pwia->Draw("sameshistE0");
   leg_Q2->AddEntry(simc_Q2_pwia,"SIMC: PWIA");
   leg_Q2->Draw();
     
   ck2->cd(2);
   data_omega->GetXaxis()->SetTitle("Energy Transfer, #omega [GeV]");
   data_omega->GetXaxis()->CenterTitle();  
   simc_omega_fsi->Draw("histE0");
   data_omega->Draw("samesE0");
   leg_om->AddEntry(data_omega,"Data", "f");
   leg_om->AddEntry(simc_omega_fsi,"SIMC: FSI");
   simc_omega_pwia->Draw("sameshistE0");
   leg_om->AddEntry(simc_omega_pwia,"SIMC: PWIA");
   leg_om->Draw();
   
   ck2->cd(3);
   data_xbj->GetXaxis()->SetTitle("BjorkenX,  X_{bj} ");
   data_xbj->GetXaxis()->CenterTitle();
   simc_xbj_fsi->Draw("histE0");
   data_xbj->Draw("samesE0");
   leg_xbj->AddEntry(data_xbj,"Data","f");
   leg_xbj->AddEntry(simc_xbj_fsi,"SIMC: FSI");
   simc_xbj_pwia->Draw("sameshistE0");
   leg_xbj->AddEntry(simc_xbj_pwia,"SIMC: PWIA");
   leg_xbj->Draw();

   ck2->cd(4);
   data_W2->GetXaxis()->SetTitle("Invariant Mass , W2 [GeV]");
   data_W2->GetXaxis()->CenterTitle();
   simc_W2_fsi->Draw("histE0");
   data_W2->Draw("samesE0");
   leg_W2->AddEntry(data_W2,"Data", "f");
   leg_W2->AddEntry(simc_W2_fsi,"SIMC: FSI");
   simc_W2_pwia->Draw("sameshistE0");
   leg_W2->AddEntry(simc_W2_pwia,"SIMC: PWIA");
   leg_W2->Draw();
   
   ck2->cd(5);
   data_th_elec->GetXaxis()->SetTitle("Electron Scatt. Angle, #theta_{e} [deg]");
   data_th_elec->GetXaxis()->CenterTitle();
   simc_th_elec_fsi->Draw("histE0");
   data_th_elec->Draw("samesE0");
   leg_the->AddEntry(data_th_elec,"Data","f");
   leg_the->AddEntry(simc_th_elec_fsi,"SIMC: FSI");
   simc_th_elec_pwia->Draw("sameshistE0");
   leg_the->AddEntry(simc_th_elec_pwia,"SIMC: PWIA");
   leg_the->Draw();
   
   ck2->cd(6);
   data_kf->GetXaxis()->SetTitle("Electron Final Momentum, k_{f} [GeV/c] ");
   data_kf->GetXaxis()->CenterTitle();   
   simc_kf_fsi->Draw("histE0");
   data_kf->Draw("samesE0");
   leg_kf->AddEntry(data_kf,"Data","f");
   leg_kf->AddEntry(simc_kf_fsi,"SIMC: FSI");
   simc_kf_pwia->Draw("sameshistE0");
   leg_kf->AddEntry(simc_kf_pwia,"SIMC: PWIA");
   leg_kf->Draw();

   ck2->cd(7);
   data_thq->GetXaxis()->SetTitle("q-vector Angle, #theta_{q} [deg]");
   data_thq->GetXaxis()->CenterTitle();
   simc_thq_fsi->Draw("histE0");
   data_thq->Draw("samesE0");
   leg_thq->AddEntry(data_thq,"Data", "f");
   leg_thq->AddEntry(simc_thq_fsi,"SIMC: FSI");
   simc_thq_pwia->Draw("sameshistE0");
   leg_thq->AddEntry(simc_thq_pwia,"SIMC: PWIA");
   leg_thq->Draw();

   ck2->cd(8);
   data_q->GetXaxis()->SetTitle("q-Vector Magnitude, |q| [GeV]");
   data_q->GetXaxis()->CenterTitle();
   simc_q_fsi->Draw("histE0");
   data_q->Draw("samesE0");
   leg_q->AddEntry(data_q,"Data", "f");
   leg_q->AddEntry(simc_q_fsi,"SIMC: FSI");
   simc_q_pwia->Draw("sameshistE0");
   leg_q->AddEntry(simc_q_pwia,"SIMC: PWIA");
   leg_q->Draw();

   ck2->SaveAs(Form("Kinematics2_pm%d_set%d.pdf", pm, set));                                                                   


   
   //Kinematics 3: Proton Kinematics
   
   auto leg_Pf = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_thp = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_Kp = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_Ep = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_Kn = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_En = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_thpq = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_thnq = new TLegend(0.1,0.8,0.28,0.9);


   //Create A Canvas to store kinematic variable comparisons
   TCanvas *ck3 = new TCanvas("ck3", "Kinematics3", 5000, 3000);
   
   ck3->Divide(4,2);
 
   ck3->cd(1);
   data_Pf->GetXaxis()->SetTitle("Proton Momentum, P_{p} [GeV]");
   data_Pf->GetXaxis()->CenterTitle();
   simc_Pf_fsi->Draw("histE0");
   data_Pf->Draw("samesE0");
   leg_Pf->AddEntry(data_Pf,"Data", "f");
   leg_Pf->AddEntry(simc_Pf_fsi,"SIMC: FSI");
   simc_Pf_pwia->Draw("sameshistE0");
   leg_Pf->AddEntry(simc_Pf_pwia,"SIMC: PWIA");
   leg_Pf->Draw();

   ck3->cd(2);
   data_th_prot->GetXaxis()->SetTitle("Proton Scatt. Angle, #theta_{p} [deg]");
   data_th_prot->GetXaxis()->CenterTitle();
   simc_th_prot_fsi->Draw("histE0");
   data_th_prot->Draw("samesE0");
   leg_thp->AddEntry(data_th_prot,"Data", "f");
   leg_thp->AddEntry(simc_th_prot_fsi,"SIMC: FSI");
   simc_th_prot_pwia->Draw("sameshistE0");
   leg_thp->AddEntry(simc_th_prot_pwia,"SIMC: PWIA");
   leg_thp->Draw();

   ck3->cd(3);
   data_Kp->GetXaxis()->SetTitle("Proton Kin. Energy, K_{p} [GeV]");
   data_Kp->GetXaxis()->CenterTitle();
   simc_Kp_fsi->Draw("histE0");
   data_Kp->Draw("samesE0");
   leg_Kp->AddEntry(data_Kp,"Data", "f");
   leg_Kp->AddEntry(simc_Kp_fsi,"SIMC: FSI");
   simc_Kp_pwia->Draw("sameshistE0");
   leg_Kp->AddEntry(simc_Kp_pwia,"SIMC: PWIA");
   leg_Kp->Draw();

   ck3->cd(4);
   data_Ep->GetXaxis()->SetTitle("Proton Energy, E_{p} [GeV]");
   data_Ep->GetXaxis()->CenterTitle();
   simc_Ep_fsi->Draw("histE0");
   data_Ep->Draw("samesE0");
   leg_Ep->AddEntry(data_Ep,"Data", "f");
   leg_Ep->AddEntry(simc_Ep_fsi,"SIMC: FSI");
   simc_Ep_pwia->Draw("sameshistE0");
   leg_Ep->AddEntry(simc_Ep_pwia,"SIMC: PWIA");
   leg_Ep->Draw();

   ck3->cd(5);
   data_Kn->GetXaxis()->SetTitle("Neutron Kin. Energy, K_{n} [GeV]");
   data_Kn->GetXaxis()->CenterTitle();
   simc_Kn_fsi->Draw("histE0");
   data_Kn->Draw("samesE0");
   leg_Kn->AddEntry(data_Kn,"Data", "f");
   leg_Kn->AddEntry(simc_Kn_fsi,"SIMC: FSI");
   simc_Kn_pwia->Draw("sameshistE0");
   leg_Kn->AddEntry(simc_Kn_pwia,"SIMC: PWIA");
   leg_Kn->Draw();

   ck3->cd(6);
   data_En->GetXaxis()->SetTitle("Neutron Energy, E_{n} [GeV]");
   data_En->GetXaxis()->CenterTitle();
   simc_En_fsi->Draw("histE0");
   data_En->Draw("samesE0");
   leg_En->AddEntry(data_En,"Data", "f");
   leg_En->AddEntry(simc_En_fsi,"SIMC: FSI");
   simc_En_pwia->Draw("sameshistE0");
   leg_En->AddEntry(simc_En_pwia,"SIMC: PWIA");
   leg_En->Draw();


   ck3->cd(7);
   data_thpq->GetXaxis()->SetTitle("(Proton, qVec.) Angle, #theta_{pq} [deg]");
   data_thpq->GetXaxis()->CenterTitle();
   simc_thpq_fsi->Draw("histE0");
   data_thpq->Draw("samesE0");
   leg_thpq->AddEntry(data_thpq,"Data", "f");
   leg_thpq->AddEntry(simc_thpq_fsi,"SIMC: FSI");
   simc_thpq_pwia->Draw("sameshistE0");
   leg_thpq->AddEntry(simc_thpq_pwia,"SIMC: PWIA");
   leg_thpq->Draw();

   ck3->cd(8);
   data_thnq->GetXaxis()->SetTitle("(Neutron, qVec.) Angle, #theta_{nq} [deg]");
   data_thnq->GetXaxis()->CenterTitle();
   simc_thnq_fsi->Draw("histE0");
   data_thnq->Draw("samesE0");
   leg_thnq->AddEntry(data_thnq,"Data", "f");
   leg_thnq->AddEntry(simc_thnq_fsi,"SIMC: FSI");
   simc_thnq_pwia->Draw("sameshistE0");
   leg_thnq->AddEntry(simc_thnq_pwia,"SIMC: PWIA");
   leg_thnq->Draw();

   ck3->SaveAs(Form("Kinematics3_pm%d_set%d.pdf", pm, set));                                                                   
                                                               

 //-----------------PLOT TARGET  Variables SIMC/Data comparison-----------------------

  //Set Legend
   auto leghxt = new TLegend(0.1,0.8,0.28,0.9);
   auto leghyt = new TLegend(0.1,0.8,0.28,0.9);
   auto leghzt = new TLegend(0.1,0.8,0.28,0.9);
   
   auto legpxt = new TLegend(0.1,0.8,0.28,0.9);
   auto legpyt = new TLegend(0.1,0.8,0.28,0.9);
   auto legpzt = new TLegend(0.1,0.8,0.28,0.9);


   TCanvas *c4a = new TCanvas("c4a", "HMS Target Variables", 5000, 3000);
   c4a->Divide(3,1);

   c4a->cd(1);
   simc_xtar_fsi->Draw("histE0");
   data_xtarH->Draw("samesE0");
   leghxt->AddEntry(data_xtarH,"Data","f");
   leghxt->AddEntry(simc_xtar_fsi,"SIMC: FSI");
   simc_xtar_pwia->Draw("sameshistE0");
   leghxt->AddEntry(simc_xtar_pwia,"SIMC: PWIA");
   leghxt->Draw();
  
   c4a->cd(2);
   simc_ytarH_fsi->Draw("histE0");
   data_ytarH->Draw("samesE0");
   leghyt->AddEntry(data_ytarH,"Data","f");
   leghyt->AddEntry(simc_ytarH_fsi,"SIMC: FSI");
   simc_ytarH_pwia->Draw("sameshistE0");
   leghyt->AddEntry(simc_ytarH_pwia,"SIMC: PWIA");
   leghyt->Draw();

   c4a->cd(3);
   simc_ztarH_fsi->Draw("histE0");
   data_ztarH->Draw("samesE0");
   leghzt->AddEntry(data_ztarH,"Data","f");
   leghzt->AddEntry(simc_ztarH_fsi,"SIMC: FSI");
   simc_ztarH_pwia->Draw("sameshistE0");
   leghzt->AddEntry(simc_ztarH_pwia,"SIMC: PWIA");
   leghzt->Draw();
  
   c4a->SaveAs(Form("hArm_TargVar_pm%d_set%d.pdf", pm, set));                                                                                              

   TCanvas *c4b = new TCanvas("c4b", "SHMS Target Variables", 5000, 3000);
   c4b->Divide(3,1);

   c4b->cd(1);
   simc_xtar_fsi->Draw("histE0");
   data_xtarP->Draw("samesE0");
   legpxt->AddEntry(data_xtarP,"Data","f");
   legpxt->AddEntry(simc_xtar_fsi,"SIMC: FSI");
   simc_xtar_pwia->Draw("sameshistE0");
   legpxt->AddEntry(simc_xtar_pwia,"SIMC: PWIA");
   legpxt->Draw();
  
   c4b->cd(2);
   simc_ytarP_fsi->Draw("histE0");
   data_ytarP->Draw("samesE0");
   legpyt->AddEntry(data_ytarP,"Data","f");
   legpyt->AddEntry(simc_ytarP_fsi,"SIMC: FSI");
   simc_ytarP_pwia->Draw("sameshistE0");
   legpyt->AddEntry(simc_ytarP_pwia,"SIMC: PWIA");
   legpyt->Draw();

   c4b->cd(3);
   simc_ztarP_fsi->Draw("histE0");
   data_ztarP->Draw("samesE0");
   legpzt->AddEntry(data_ztarP,"Data","f");
   legpzt->AddEntry(simc_ztarP_fsi,"SIMC: FSI");
   simc_ztarP_pwia->Draw("sameshistE0");
   legpzt->AddEntry(simc_ztarP_pwia,"SIMC: PWIA");
   legpzt->Draw();
  
   c4b->SaveAs(Form("pArm_TargVar_pm%d_set%d.pdf", pm, set));      
   //--------PLOT HADRON ARM QUANTITIES--------


   
   //-----------------PLOT Target Reconstructed Variables SIMC/Data comparison-----------------------
 
   //Set Legend
   auto htr_l1 = new TLegend(0.1,0.8,0.28,0.9);
   auto htr_l2 = new TLegend(0.1,0.8,0.28,0.9);
   auto htr_l3 = new TLegend(0.1,0.8,0.28,0.9);
   auto htr_l4 = new TLegend(0.1,0.8,0.28,0.9);
   
   //Create A Canvas to store Target Recon. variable comparisons in HADRON ARM
   
   TCanvas *htr = new TCanvas("htr", "Hadron Arm: Target Reconstruction", 5000, 3000);
   htr->Divide(2,2);

   htr->cd(1);
   simc_hytar_fsi->Draw("histE0");
   data_hytar->Draw("samesE0");
   htr_l1->AddEntry(data_hytar,"Data","f");
   htr_l1->AddEntry(simc_hytar_fsi,"SIMC: FSI");
   simc_hytar_pwia->Draw("sameshistE0");
   htr_l1->AddEntry(simc_hytar_pwia,"SIMC: PWIA");
   htr_l1->Draw();

   htr->cd(2);
   simc_hxptar_fsi->Draw("histE0");
   data_hxptar->Draw("samesE0");
   htr_l2->AddEntry(data_hxptar,"Data", "f");
   htr_l2->AddEntry(simc_hxptar_fsi,"SIMC: FSI");
   simc_hxptar_pwia->Draw("sameshistE0");
   htr_l2->AddEntry(simc_hxptar_pwia,"SIMC: PWIA");
   htr_l2->Draw();

   htr->cd(3);
   simc_hyptar_fsi->Draw("histE0");
   data_hyptar->Draw("samesE0");
   htr_l3->AddEntry(data_hyptar,"Data", "f");
   htr_l3->AddEntry(simc_hyptar_fsi,"SIMC: FSI");
   simc_hyptar_pwia->Draw("sameshistE0");
   htr_l3->AddEntry(simc_hyptar_pwia,"SIMC: PWIA");
   htr_l3->Draw();
     
   htr->cd(4);
   simc_hdelta_fsi->Draw("histE0");
   data_hdelta->Draw("samesE0");
   htr_l4->AddEntry(data_hdelta,"Data", "f");
   htr_l4->AddEntry(simc_hdelta_fsi,"SIMC: FSI");
   simc_hdelta_pwia->Draw("sameshistE0");
   htr_l4->AddEntry(simc_hdelta_pwia,"SIMC: PWIA");
   htr_l4->Draw();

   htr->SaveAs(Form("hArm_TargRecon_pm%d_set%d.pdf", pm, set));

   //------------------------------------------------------------------------------

   
   //-----------------PLOT FOCAL PLANE  Variables SIMC/Data comparison-----------------------

   //Set Legend
   auto hfp_l1 = new TLegend(0.1,0.8,0.28,0.9);
   auto hfp_l2 = new TLegend(0.1,0.8,0.28,0.9);
   auto hfp_l3 = new TLegend(0.1,0.8,0.28,0.9);
   auto hfp_l4 = new TLegend(0.1,0.8,0.28,0.9);

   TCanvas *hfp = new TCanvas("hfp", "Hadron Arm: Focal Plane", 5000, 3000);
   hfp->Divide(2,2);

   hfp->cd(1);
   simc_hxfp_fsi->Draw("histE0");
   data_hxfp->Draw("samesE0");
   hfp_l1->AddEntry(data_hxfp,"Data","f");
   hfp_l1->AddEntry(simc_hxfp_fsi,"SIMC: FSI");
   simc_hxfp_pwia->Draw("sameshistE0");
   hfp_l1->AddEntry(simc_hxfp_pwia,"SIMC: PWIA");
   hfp_l1->Draw();
   
   hfp->cd(2);
   simc_hyfp_fsi->Draw("histE0");
   data_hyfp->Draw("samesE0");
   hfp_l2->AddEntry(data_hyfp,"Data", "f");
   hfp_l2->AddEntry(simc_hyfp_fsi,"SIMC: FSI");
   simc_hyfp_pwia->Draw("sameshistE0");
   hfp_l2->AddEntry(simc_hyfp_pwia,"SIMC: PWIA");
   hfp_l2->Draw();

   hfp->cd(3);
   simc_hxpfp_fsi->Draw("histE0");
   data_hxpfp->Draw("samesE0");
   hfp_l3->AddEntry(data_hxpfp,"Data", "f");
   hfp_l3->AddEntry(simc_hxpfp_fsi,"SIMC: FSI");
   simc_hxpfp_pwia->Draw("sameshistE0");
   hfp_l3->AddEntry(simc_hxpfp_pwia,"SIMC: PWIA");
   hfp_l3->Draw();
     
   hfp->cd(4);
   simc_hypfp_fsi->Draw("histE0");
   data_hypfp->Draw("samesE0");
   hfp_l4->AddEntry(data_hypfp,"Data", "f");
   hfp_l4->AddEntry(simc_hypfp_fsi,"SIMC: FSI");
   simc_hypfp_pwia->Draw("sameshistE0");
   hfp_l4->AddEntry(simc_hypfp_pwia,"SIMC: PWIA");
   hfp_l4->Draw();

   hfp->SaveAs(Form("hArm_FocalPlane_pm%d_set%d.pdf", pm, set));                                                                                   

   //----------------------------------------------------------- 
 
     //Set Legend
   auto leg_pmr = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_thnqr = new TLegend(0.1,0.8,0.28,0.9);

  
   //-----------------PLOT Target Reconstructed Variables SIMC/Data comparison-----------------------

   //Create A Canvas to store Target Recon. variable comparisons in HADRON ARM
   
   TCanvas *cR = new TCanvas("cR", "Data/SIMC Ratio", 5000, 3000);
   cR->Divide(1,2);

   cR->cd(1);
   dataPm_fsi->Divide(simc_Pm_fsi);
   dataPm_fsi->SetLineColor(kRed);
   dataPm_fsi->SetLineWidth(2);
   dataPm_fsi->SetMarkerStyle(21);
   dataPm_fsi->SetMarkerColor(kRed);
   dataPm_fsi->Draw();
   dataPm_pwia->Divide(simc_Pm_pwia);
   dataPm_pwia->SetLineColor(kBlue);
   dataPm_pwia->SetLineWidth(2);
   dataPm_pwia->SetMarkerStyle(21);
   dataPm_pwia->SetMarkerColor(kBlue);
   dataPm_pwia->Draw("sames");
   dataPm_fsi->GetXaxis()->SetRangeUser(0.3,1.3);
   dataPm_fsi->GetYaxis()->SetRangeUser(0.3,1.7);
   leg_pmr->AddEntry(dataPm_fsi,"FSI Ratio","f");
   leg_pmr->AddEntry(dataPm_pwia,"PWIA Ratio","f");
   leg_pmr->Draw();

   cR->cd(2);
   data_thnq_fsi->Divide(simc_thnq_fsi);
   data_thnq_fsi->SetLineColor(kRed);
   data_thnq_fsi->SetLineWidth(2);
   data_thnq_fsi->SetMarkerStyle(21);
   data_thnq_fsi->SetMarkerColor(kRed);
   data_thnq_fsi->Draw();
   data_thnq_pwia->Divide(simc_thnq_pwia);
   data_thnq_pwia->SetLineColor(kBlue);
   data_thnq_pwia->SetLineWidth(2);
   data_thnq_pwia->SetMarkerStyle(21);
   data_thnq_pwia->SetMarkerColor(kBlue);
   data_thnq_pwia->Draw("sames");
   data_thnq_fsi->GetXaxis()->SetRangeUser(0,100);
   data_thnq_fsi->GetYaxis()->SetRangeUser(0.3,1.7);
   leg_thnqr->AddEntry(data_thnq_fsi,"FSI Ratio","f");
   leg_thnqr->AddEntry(data_thnq_pwia,"PWIA Ratio","f");
   leg_thnqr->Draw();


   cR->SaveAs(Form("YieldRatio_pm%d_set%d.pdf", pm, set));                                        



}
