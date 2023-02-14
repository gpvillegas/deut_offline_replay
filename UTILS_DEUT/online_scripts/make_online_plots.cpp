// Script to either plot ONLY DATA or make comparison between DATA and SIMC 
// Histogram objects are retrieved from pre-existing ROOTfiles with pre-determined
// histogram names

/* histogram categories:
1) HMS/SHMS focal plane,
2) HMS/SHMS Reconstructed,
3) Kinematics-1,2  
4) Target Vertex
*/

void make_online_plots(int run=0, int evt=0, Bool_t simc_exist=0, TString tgt_type="", TString ana_type="", TString ana_cut="", TString data_file_path="", TString simc_file_path="", Bool_t draw_norm=1)
{

  gROOT->SetBatch(kTRUE);  
  gStyle->SetOptStat(1001111);
  
  //Pre-defined SIMC/data root file names containing histogram object to comapare

  //TString simc_filename =  Form("../heep_simc_histos_%d_rad.root", run);                      
  //TString data_filename = Form("../heep_data_histos_%d_combined.root",run); 

  TString outPDF=Form("DEUT_OUTPUT/PDF/deut_output_%s_%d_%d.pdf", ana_type.Data(), run, evt);
  
  Bool_t data_exist = !gSystem->AccessPathName( data_file_path.Data() );
  if(!data_exist){
    cout << Form("data file: %s does NOT exist. Exit. ", data_file_path.Data() ) << endl;
    gSystem->Exit(0);
  }
  
  cout << "simc_exist ? " << simc_exist << endl;
  //Where to store plots
  string plots_dir = "./";
  string plots_path;
  
  //Open SIMC/data ROOT files;
  TFile *simc_file = NULL;
  TFile *data_file = NULL;
  if(simc_exist) simc_file = new TFile(simc_file_path.Data());
  data_file = new TFile(data_file_path.Data());

 
  //---------------Target ----------------
  //Define SIMC histos ('h'-->hadron arm,  'e'-->electron arm)
  
  TH1F *simc_xtarH =  0;
  TH1F *simc_ytarH =  0;
  TH1F *simc_ztarH =  0;

  TH1F *simc_xtarP =  0;
  TH1F *simc_ytarP =  0;
  TH1F *simc_ztarP =  0;  

  //Define data histos
  TH1F *data_xtarH = 0;
  TH1F *data_ytarH = 0;
  TH1F *data_ztarH = 0;

  TH1F *data_xtarP = 0;                                                                                                     
  TH1F *data_ytarP = 0;                                                                                                                                   
  TH1F *data_ztarP = 0; 

  //---------------Target Reconstruction Variables----------------
  //Define SIMC histos ('h'-->hadron arm,  'e'-->electron arm)
  TH1F *simc_eytar =  0;
  TH1F *simc_exptar =  0;
  TH1F *simc_eyptar =  0;
  TH1F *simc_edelta =  0;

  TH1F *simc_hytar =  0;
  TH1F *simc_hxptar =  0;
  TH1F *simc_hyptar =  0;
  TH1F *simc_hdelta =  0;

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
  TH1F *simc_exfp =  0;
  TH1F *simc_eyfp =  0;
  TH1F *simc_expfp =  0;
  TH1F *simc_eypfp =  0;

  TH1F *simc_hxfp =  0;
  TH1F *simc_hyfp =  0;
  TH1F *simc_hxpfp =  0;
  TH1F *simc_hypfp =  0;
  
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
  TH1F *simc_Q2 =  0;
  TH1F *simc_nu =  0;
  TH1F *simc_W =  0;
  TH1F *simc_thq = 0;

  TH1F *simc_xbj = 0;
  TH1F *simc_the = 0;                                  
  TH1F *simc_kf = 0;  
  TH1F *simc_Em = 0;
  TH1F *simc_MM = 0;

  //Kinematics 2
  TH1F *simc_Pm = 0;
  TH1F *simc_Pf = 0;
  TH1F *simc_thx = 0;
  TH1F *simc_q = 0;    //q-vector magnitude
  TH1F *simc_thxq = 0;
  TH1F *simc_thrq = 0;
  TH1F *simc_Pmx = 0;
  TH1F *simc_Pmy = 0;
  TH1F *simc_Pmz = 0;

  // 2d kinematics
  TH2F * simc_Em_nuc_vs_Pm = 0;
  TH2F * simc_Em_src_vs_Pm = 0;
  
  //Define data histos
  TH1F *data_Q2 =  0;
  TH1F *data_nu =  0;
  TH1F *data_W =  0;
  TH1F *data_thq = 0;

  TH1F *data_xbj = 0;
  TH1F *data_the = 0;
  TH1F *data_kf = 0;
  TH1F *data_Em = 0;
  TH1F *data_MM = 0;

   //Kinematics 2
  TH1F *data_Pm = 0;
  TH1F *data_Pf = 0;
  TH1F *data_thx = 0;
  TH1F *data_q = 0;    //q-vector magnitude
  TH1F *data_thxq = 0;
  TH1F *data_thrq = 0;
  TH1F *data_Pmx = 0;
  TH1F *data_Pmy = 0;
  TH1F *data_Pmz = 0;
  // 2d kinematics
  TH2F * data_Em_nuc_vs_Pm = 0;
  TH2F * data_Em_src_vs_Pm = 0;
  
  //---------------------------------------------------------------

  
  //---------------- SELECTED DATA: TOTAL = SIGNAL + BACKGROUND PLOTS -----------------

  //  electron-proton Coincidence Time
  TH1F *data_ep_ctime_total = 0;  // full coin. time spectrum
  TH1F *data_ep_ctime_real = 0;  // main coin. peak selected
  TH1F *data_ep_ctime_rand = 0;  // out-of-time ("accidental") randoms selected 

  // Invariant Mass
  TH1F *data_W_total = 0;   
  TH1F *data_W_real = 0;  // random coincidence subtracted
  TH1F *data_W_rand = 0;  // random coincidences selected

  // Missing Mass
  TH1F *data_MM_total = 0;   
  TH1F *data_MM_real = 0;  // random coincidence subtracted
  TH1F *data_MM_rand = 0;  // random coincidences selected

  // Missing Momentum
  TH1F *data_Pm_total = 0;   
  TH1F *data_Pm_real = 0;  // random coincidence subtracted
  TH1F *data_Pm_rand = 0;  // random coincidences selected
  
  // Missing Energy
  TH1F *data_Em_total = 0;   
  TH1F *data_Em_real = 0;  // random coincidence subtracted
  TH1F *data_Em_rand = 0;  // random coincidences selected
  
  
  //-----------------------------------------------------------------------------------

  
  if(simc_exist) {

    
    //change to simc_file
    simc_file->cd();
    
    
    //----------Get Target Histograms------------------
    //Get Histogram objects from SIMC rootfile
    simc_file->GetObject("accp_plots/H_htar_x", simc_xtarH);
    simc_file->GetObject("accp_plots/H_htar_y", simc_ytarH);
    simc_file->GetObject("accp_plots/H_htar_z", simc_ztarH);
    
    simc_file->GetObject("accp_plots/H_etar_x", simc_xtarP);
    simc_file->GetObject("accp_plots/H_etar_y", simc_ytarP);  
    simc_file->GetObject("accp_plots/H_etar_z", simc_ztarP); 
    
    //Set SIMC Histo Aesthetics
    simc_xtarH->SetLineColor(kRed);
    simc_xtarH->SetLineWidth(2);
    
    simc_xtarP->SetLineColor(kRed);
    simc_xtarP->SetLineWidth(2);

    simc_ytarH->SetLineColor(kRed);
    simc_ytarH->SetLineWidth(2);
    simc_ztarH->SetLineColor(kRed);
    simc_ztarH->SetLineWidth(2);
    
    simc_ytarP->SetLineColor(kRed);          
    simc_ytarP->SetLineWidth(2);           
    simc_ztarP->SetLineColor(kRed);                        
    simc_ztarP->SetLineWidth(2);  
    
  }
  
  //change to data_file
  data_file->cd();

  //Get Histogram objects from data rootfile
  data_file->GetObject("accp_plots/H_htar_x", data_xtarH);
  data_file->GetObject("accp_plots/H_htar_y", data_ytarH);
  data_file->GetObject("accp_plots/H_htar_z", data_ztarH);

  data_file->GetObject("accp_plots/H_etar_x", data_xtarP); 
  data_file->GetObject("accp_plots/H_etar_y", data_ytarP);                        
  data_file->GetObject("accp_plots/H_etar_z", data_ztarP); 
  
  //Set data Histo Aesthetics
  data_xtarH->SetFillColorAlpha(kBlue, 0.35);
  data_xtarH->SetFillStyle(3004);
  data_ytarH->SetFillColorAlpha(kBlue, 0.35);
  data_ytarH->SetFillStyle(3004);
  data_ztarH->SetFillColorAlpha(kBlue, 0.35);
  data_ztarH->SetFillStyle(3004);

  data_xtarP->SetFillColorAlpha(kBlue, 0.35);         
  data_xtarP->SetFillStyle(3004); 
  data_ytarP->SetFillColorAlpha(kBlue, 0.35);                                  
  data_ytarP->SetFillStyle(3004);                               
  data_ztarP->SetFillColorAlpha(kBlue, 0.35);             
  data_ztarP->SetFillStyle(3004);  

  //-----------------------------------------------------------------


  //---------------------------------------------------------------

  if(simc_exist) {
    //change to simc_file
    simc_file->cd();
    
    //----------Get Target Reconstructed Histograms------------------
    //Get Histogram objects from SIMC rootfile
    simc_file->GetObject("accp_plots/H_eytar", simc_eytar);
    simc_file->GetObject("accp_plots/H_exptar", simc_exptar);
    simc_file->GetObject("accp_plots/H_eyptar", simc_eyptar);
    simc_file->GetObject("accp_plots/H_edelta", simc_edelta);
    
    simc_file->GetObject("accp_plots/H_hytar", simc_hytar);
    simc_file->GetObject("accp_plots/H_hxptar", simc_hxptar);
    simc_file->GetObject("accp_plots/H_hyptar", simc_hyptar);
    simc_file->GetObject("accp_plots/H_hdelta", simc_hdelta);

    
    //Set SIMC Histo Aesthetics
    simc_eytar->SetLineColor(kRed);
    simc_eytar->SetLineWidth(2);
    simc_exptar->SetLineColor(kRed);
    simc_exptar->SetLineWidth(2);
    simc_eyptar->SetLineColor(kRed);
    simc_eyptar->SetLineWidth(2);
    simc_edelta->SetLineColor(kRed);
    simc_edelta->SetLineWidth(2);
    
    simc_hytar->SetLineColor(kRed);
    simc_hytar->SetLineWidth(2);
    simc_hxptar->SetLineColor(kRed);
    simc_hxptar->SetLineWidth(2);
    simc_hyptar->SetLineColor(kRed);
    simc_hyptar->SetLineWidth(2);
    simc_hdelta->SetLineColor(kRed);
    simc_hdelta->SetLineWidth(2);
    
  }
  
  //change to data_file
  data_file->cd();

  //Get Histogram objects from data rootfile
  data_file->GetObject("accp_plots/H_eytar", data_eytar);
  data_file->GetObject("accp_plots/H_exptar", data_exptar);
  data_file->GetObject("accp_plots/H_eyptar", data_eyptar);
  data_file->GetObject("accp_plots/H_edelta", data_edelta);
  
  data_file->GetObject("accp_plots/H_hytar", data_hytar);
  data_file->GetObject("accp_plots/H_hxptar", data_hxptar);
  data_file->GetObject("accp_plots/H_hyptar", data_hyptar);
  data_file->GetObject("accp_plots/H_hdelta", data_hdelta);

  //Set data Histo Aesthetics
  data_eytar->SetFillColorAlpha(kBlue, 0.35);
  data_eytar->SetFillStyle(3004);
  data_exptar->SetFillColorAlpha(kBlue, 0.35);
  data_exptar->SetFillStyle(3004);
  data_eyptar->SetFillColorAlpha(kBlue, 0.35);
  data_eyptar->SetFillStyle(3004);
  data_edelta->SetFillColorAlpha(kBlue, 0.35);
  data_edelta->SetFillStyle(3004);

  data_hytar->SetFillColorAlpha(kBlue, 0.35);
  data_hytar->SetFillStyle(3004);
  data_hxptar->SetFillColorAlpha(kBlue, 0.35);
  data_hxptar->SetFillStyle(3004);
  data_hyptar->SetFillColorAlpha(kBlue, 0.35);
  data_hyptar->SetFillStyle(3004);
  data_hdelta->SetFillColorAlpha(kBlue, 0.35);
  data_hdelta->SetFillStyle(3004);

  //-----------------------------------------------------------------

  

  //---------------Get FOCAL PLANE Histograms------------------------
  if(simc_exist) {
   //change to simc_file
  simc_file->cd();

  //Get Histogram objects from SIMC rootfile
  simc_file->GetObject("accp_plots/H_exfp", simc_exfp);
  simc_file->GetObject("accp_plots/H_eyfp", simc_eyfp);
  simc_file->GetObject("accp_plots/H_expfp", simc_expfp);
  simc_file->GetObject("accp_plots/H_eypfp", simc_eypfp);

  simc_file->GetObject("accp_plots/H_hxfp", simc_hxfp);
  simc_file->GetObject("accp_plots/H_hyfp", simc_hyfp);
  simc_file->GetObject("accp_plots/H_hxpfp", simc_hxpfp);
  simc_file->GetObject("accp_plots/H_hypfp", simc_hypfp);
  //Set SIMC Histo Aesthetics
  simc_exfp->SetLineColor(kRed);
  simc_exfp->SetLineWidth(2);
  simc_eyfp->SetLineColor(kRed);
  simc_eyfp->SetLineWidth(2);
  simc_expfp->SetLineColor(kRed);
  simc_expfp->SetLineWidth(2);
  simc_eypfp->SetLineColor(kRed);
  simc_eypfp->SetLineWidth(2);
  
  simc_hxfp->SetLineColor(kRed);
  simc_hxfp->SetLineWidth(2);
  simc_hyfp->SetLineColor(kRed);
  simc_hyfp->SetLineWidth(2);
  simc_hxpfp->SetLineColor(kRed);
  simc_hxpfp->SetLineWidth(2);
  simc_hypfp->SetLineColor(kRed);
  simc_hypfp->SetLineWidth(2);
  }
  
  //change to data_file
  data_file->cd();

  //Get Histogram objects from data rootfile
  data_file->GetObject("accp_plots/H_exfp", data_exfp);
  data_file->GetObject("accp_plots/H_eyfp", data_eyfp);
  data_file->GetObject("accp_plots/H_expfp", data_expfp);
  data_file->GetObject("accp_plots/H_eypfp", data_eypfp);

  data_file->GetObject("accp_plots/H_hxfp", data_hxfp);
  data_file->GetObject("accp_plots/H_hyfp", data_hyfp);
  data_file->GetObject("accp_plots/H_hxpfp", data_hxpfp);
  data_file->GetObject("accp_plots/H_hypfp", data_hypfp);
  //Set data Histo Aesthetics
  data_exfp->SetFillColorAlpha(kBlue, 0.35);
  data_exfp->SetFillStyle(3004);
  data_eyfp->SetFillColorAlpha(kBlue, 0.35);
  data_eyfp->SetFillStyle(3004);
  data_expfp->SetFillColorAlpha(kBlue, 0.35);
  data_expfp->SetFillStyle(3004);
  data_eypfp->SetFillColorAlpha(kBlue, 0.35);
  data_eypfp->SetFillStyle(3004);

  data_hxfp->SetFillColorAlpha(kBlue, 0.35);
  data_hxfp->SetFillStyle(3004);
  data_hyfp->SetFillColorAlpha(kBlue, 0.35);
  data_hyfp->SetFillStyle(3004);
  data_hxpfp->SetFillColorAlpha(kBlue, 0.35);
  data_hxpfp->SetFillStyle(3004);
  data_hypfp->SetFillColorAlpha(kBlue, 0.35);
  data_hypfp->SetFillStyle(3004);

  //--------------------------------------------------------------
  
  //------------------Get KINEMATICS VARIABLES--------------------

  if(simc_exist) {
   //change to simc_file
  simc_file->cd();

  //Get Histogram objects from SIMC rootfile
  simc_file->GetObject("kin_plots/H_Q2", simc_Q2);
  simc_file->GetObject("kin_plots/H_nu", simc_nu);
  simc_file->GetObject("kin_plots/H_W", simc_W);
  simc_file->GetObject("kin_plots/H_thq", simc_thq);

  simc_file->GetObject("kin_plots/H_xbj", simc_xbj);
  simc_file->GetObject("kin_plots/H_the", simc_the);
  simc_file->GetObject("kin_plots/H_kf", simc_kf);
  simc_file->GetObject("kin_plots/H_Em", simc_Em);
  simc_file->GetObject("kin_plots/H_MM", simc_MM);

  simc_file->GetObject("kin_plots/H_Pm", simc_Pm);
  simc_file->GetObject("kin_plots/H_Pf", simc_Pf);
  simc_file->GetObject("kin_plots/H_thx", simc_thx);
  simc_file->GetObject("kin_plots/H_q", simc_q);
  simc_file->GetObject("kin_plots/H_thxq", simc_thxq);
  simc_file->GetObject("kin_plots/H_thrq", simc_thrq);

  simc_file->GetObject("kin_plots/H_Pmx_Lab", simc_Pmx);
  simc_file->GetObject("kin_plots/H_Pmy_Lab", simc_Pmy);
  simc_file->GetObject("kin_plots/H_Pmz_Lab", simc_Pmz);
  
  //Set SIMC Histo Aesthetics
  simc_Q2->SetLineColor(kRed);
  simc_Q2->SetLineWidth(2);
  simc_nu->SetLineColor(kRed);
  simc_nu->SetLineWidth(2);
  simc_W->SetLineColor(kRed);
  simc_W->SetLineWidth(2);
  simc_thq->SetLineColor(kRed);
  simc_thq->SetLineWidth(2);
  
  simc_xbj->SetLineColor(kRed);
  simc_xbj->SetLineWidth(2);
  simc_the->SetLineColor(kRed);
  simc_the->SetLineWidth(2);
  simc_kf->SetLineColor(kRed);
  simc_kf->SetLineWidth(2);
  simc_Em->SetLineColor(kRed);
  simc_Em->SetLineWidth(2);
  simc_MM->SetLineColor(kRed);
  simc_MM->SetLineWidth(2);
  
  simc_Pm->SetLineColor(kRed);
  simc_Pm->SetLineWidth(2);
  simc_Pf->SetLineColor(kRed);
  simc_Pf->SetLineWidth(2);
  simc_thx->SetLineColor(kRed);
  simc_thx->SetLineWidth(2);
  simc_q->SetLineColor(kRed);
  simc_q->SetLineWidth(2);
  simc_thxq->SetLineColor(kRed);
  simc_thxq->SetLineWidth(2);
  simc_thrq->SetLineColor(kRed);
  simc_thrq->SetLineWidth(2);
  
  simc_Pmx->SetLineColor(kRed);
  simc_Pmx->SetLineWidth(2);
  simc_Pmy->SetLineColor(kRed);
  simc_Pmy->SetLineWidth(2);
  simc_Pmz->SetLineColor(kRed);
  simc_Pmz->SetLineWidth(2);
    
  }
  
  //change to data_file
  data_file->cd();
  
  //Get Histogram objects from data rootfile
  data_file->GetObject("kin_plots/H_Q2", data_Q2);
  data_file->GetObject("kin_plots/H_nu", data_nu);
  data_file->GetObject("kin_plots/H_W", data_W);
  data_file->GetObject("kin_plots/H_thq", data_thq);

  
  data_file->GetObject("kin_plots/H_xbj", data_xbj);
  data_file->GetObject("kin_plots/H_the", data_the);
  data_file->GetObject("kin_plots/H_kf", data_kf);

  if(tgt_type=="LH2"){
    data_file->GetObject("kin_plots/H_Em", data_Em);
  }
  else if (tgt_type!="LH2"){
    data_file->GetObject("kin_plots/H_Em_nuc", data_Em);

    data_file->GetObject("kin_plots/H_Em_nuc_vs_Pm", data_Em_nuc_vs_Pm);
    data_file->GetObject("kin_plots/H_Em_src_vs_Pm", data_Em_src_vs_Pm);
    
  }
  
  data_file->GetObject("kin_plots/H_MM", data_MM);
  data_file->GetObject("kin_plots/H_Pm", data_Pm);
  data_file->GetObject("kin_plots/H_Pf", data_Pf);
  data_file->GetObject("kin_plots/H_thx", data_thx);
  data_file->GetObject("kin_plots/H_q", data_q);
  data_file->GetObject("kin_plots/H_thxq", data_thxq);
  data_file->GetObject("kin_plots/H_thrq", data_thrq);

  data_file->GetObject("kin_plots/H_Pmx_Lab", data_Pmx);
  data_file->GetObject("kin_plots/H_Pmy_Lab", data_Pmy);
  data_file->GetObject("kin_plots/H_Pmz_Lab", data_Pmz);

  //Set data Histo Aesthetics
  data_Q2->SetFillColorAlpha(kBlue, 0.35);
  data_Q2->SetFillStyle(3004);
  data_nu->SetFillColorAlpha(kBlue, 0.35);
  data_nu->SetFillStyle(3004);
  data_W->SetFillColorAlpha(kBlue, 0.35);
  data_W->SetFillStyle(3004);
  data_thq->SetFillColorAlpha(kBlue, 0.35);
  data_thq->SetFillStyle(3004);

  data_xbj->SetFillColorAlpha(kBlue, 0.35);
  data_xbj->SetFillStyle(3004);
  data_the->SetFillColorAlpha(kBlue, 0.35);
  data_the->SetFillStyle(3004);
  data_kf->SetFillColorAlpha(kBlue, 0.35);
  data_kf->SetFillStyle(3004);
  data_Em->SetFillColorAlpha(kBlue, 0.35);
  data_Em->SetFillStyle(3004);
  data_MM->SetFillColorAlpha(kBlue, 0.35);
  data_MM->SetFillStyle(3004);
  
  data_Pm->SetFillColorAlpha(kBlue,0.35);
  data_Pm->SetFillStyle(3004);
  data_Pf->SetFillColorAlpha(kBlue,0.35);
  data_Pf->SetFillStyle(3004);
  data_thx->SetFillColorAlpha(kBlue,0.35);
  data_thx->SetFillStyle(3004);
  data_q->SetFillColorAlpha(kBlue,0.35);
  data_q->SetFillStyle(3004);
  data_thxq->SetFillColorAlpha(kBlue,0.35);
  data_thxq->SetFillStyle(3004);
  data_thrq->SetFillColorAlpha(kBlue,0.35);
  data_thrq->SetFillStyle(3004);
  
  data_Pmx->SetFillColorAlpha(kBlue,0.35);
  data_Pmx->SetFillStyle(3004);
  data_Pmy->SetFillColorAlpha(kBlue,0.35);
  data_Pmy->SetFillStyle(3004);
  data_Pmz->SetFillColorAlpha(kBlue,0.35);
  data_Pmz->SetFillStyle(3004);


  //---------------- GET SELECTED DATA: TOTAL = SIGNAL + BACKGROUND PLOTS -----------------

  //change to data_file
  data_file->cd();
  
  //Get Histogram objects from data rootfile
  data_file->GetObject("pid_plots/H_ep_ctime_total", data_ep_ctime_total);
  data_file->GetObject("pid_plots/H_ep_ctime", data_ep_ctime_real);
  data_file->GetObject("rand_plots/H_ep_ctime_rand", data_ep_ctime_rand);

  data_file->GetObject("kin_plots/H_W", data_W_total);
  data_file->GetObject("randSub_plots/H_W_rand_sub", data_W_real);
  data_file->GetObject("rand_plots/H_W_rand", data_W_rand);

  data_file->GetObject("kin_plots/H_MM", data_MM_total);
  data_file->GetObject("randSub_plots/H_MM_rand_sub", data_MM_real);
  data_file->GetObject("rand_plots/H_MM_rand", data_MM_rand);

  data_file->GetObject("kin_plots/H_Pm", data_Pm_total);
  data_file->GetObject("randSub_plots/H_Pm_rand_sub", data_Pm_real);
  data_file->GetObject("rand_plots/H_Pm_rand", data_Pm_rand);

  if(tgt_type=="LH2"){
    data_file->GetObject("kin_plots/H_Em", data_Em_total);
    data_file->GetObject("randSub_plots/H_Em_rand_sub", data_Em_real);
    data_file->GetObject("rand_plots/H_Em_rand", data_Em_rand);
  }
  else if(tgt_type!="LH2"){
    data_file->GetObject("kin_plots/H_Em_nuc", data_Em_total);
    data_file->GetObject("randSub_plots/H_Em_nuc_rand_sub", data_Em_real);
    data_file->GetObject("rand_plots/H_Em_nuc_rand", data_Em_rand);
  }
    
  //Set data Histo Aesthetics

  // coincidence time
  data_ep_ctime_total->SetFillColorAlpha(kBlue, 0.35);
  data_ep_ctime_total->SetFillStyle(3004);
  data_ep_ctime_total->SetLineColor(kBlue+2);

  data_ep_ctime_rand->SetFillColorAlpha(kGreen, 0.35);
  data_ep_ctime_rand->SetFillStyle(3005);
  data_ep_ctime_rand->SetLineColor(kGreen);
  
  data_ep_ctime_real->SetFillColorAlpha(kMagenta, 0.35);
  data_ep_ctime_real->SetFillStyle(3006);
  data_ep_ctime_real->SetLineColor(kMagenta);

  // invariant mass, W
  data_W_total->SetFillColorAlpha(kBlue, 0.35);
  data_W_total->SetFillStyle(3004);
  data_W_total->SetLineColor(kBlue+2);

  data_W_rand->SetFillColorAlpha(kGreen, 0.35);
  data_W_rand->SetFillStyle(3005);
  data_W_rand->SetLineColor(kGreen);
  
  data_W_real->SetFillColorAlpha(kMagenta, 0.35);
  data_W_real->SetFillStyle(3006);
  data_W_real->SetLineColor(kMagenta);

  // missing mass, MM
  data_MM_total->SetFillColorAlpha(kBlue, 0.35);
  data_MM_total->SetFillStyle(3004);
  data_MM_total->SetLineColor(kBlue+2);

  data_MM_rand->SetFillColorAlpha(kGreen, 0.35);
  data_MM_rand->SetFillStyle(3005);
  data_MM_rand->SetLineColor(kGreen);
  
  data_MM_real->SetFillColorAlpha(kMagenta, 0.35);
  data_MM_real->SetFillStyle(3006);
  data_MM_real->SetLineColor(kMagenta);
  
  // missing momentum
  data_Pm_total->SetFillColorAlpha(kBlue, 0.35);
  data_Pm_total->SetFillStyle(3004);
  data_Pm_total->SetLineColor(kBlue+2);

  data_Pm_rand->SetFillColorAlpha(kGreen, 0.35);
  data_Pm_rand->SetFillStyle(3005);
  data_Pm_rand->SetLineColor(kGreen);
  
  data_Pm_real->SetFillColorAlpha(kMagenta, 0.35);
  data_Pm_real->SetFillStyle(3006);
  data_Pm_real->SetLineColor(kMagenta);

  // missing energy
  data_Em_total->SetFillColorAlpha(kBlue, 0.35);
  data_Em_total->SetFillStyle(3004);
  data_Em_total->SetLineColor(kBlue+2);

  data_Em_rand->SetFillColorAlpha(kGreen, 0.35);
  data_Em_rand->SetFillStyle(3005);
  data_Em_rand->SetLineColor(kGreen);
  
  data_Em_real->SetFillColorAlpha(kMagenta, 0.35);
  data_Em_real->SetFillStyle(3006);
  data_Em_real->SetLineColor(kMagenta);

  //-----------------------------------------------------------------------------------

  gStyle->SetOptStat(0);
    
  // Create canvas to store multi-page .pdf plots
  TCanvas *c1 = new TCanvas("c1", "deut_output", 2000, 1000); 
  c1->Print(Form("deut_output_%s_%d.pdf[", ana_type.Data(), run));
  c1->Clear();
  
  //---------------- SELECTED DATA: TOTAL = SIGNAL + BACKGROUND PLOTS -----------------
  double total, total_err;
  double reals, reals_err;
  double rands, rands_err;
  
  
  double nbins;
  
  auto hctime_leg = new TLegend(0.15,0.6,0.33,0.8);
  auto hW_leg     = new TLegend(0.15,0.6,0.33,0.8);
  auto hMM_leg    = new TLegend(0.15,0.6,0.33,0.8);
  auto hPm_leg    = new TLegend(0.63,0.6,0.8,0.8);
  auto hEm_leg    = new TLegend(0.63,0.6,0.8,0.8);
  
  
  // ------- COINCIDENCE TIME -----
  if(ana_cut!="heep_singles"){
  c1->cd();
  c1->SetLogy();
  nbins = data_ep_ctime_total->GetNbinsX();  //Get total number of bins (excluding overflow) (same for total, reals randoms of same histo)

  data_ep_ctime_total->GetYaxis()->SetRangeUser(0.5, data_ep_ctime_total->GetMaximum()+1e5);
  data_ep_ctime_total->Draw("histE0");   
  data_ep_ctime_real->Draw("sameshistE0");   
  data_ep_ctime_rand->Draw("sameshistE0");   
  
  total = data_ep_ctime_total->IntegralAndError(1, nbins, total_err);
  reals = data_ep_ctime_real->IntegralAndError(1, nbins, reals_err);
  rands = data_ep_ctime_rand->IntegralAndError(1, nbins, rands_err);
  
  hctime_leg->AddEntry(data_ep_ctime_total,Form("Total   : %.3f", total),"f");
  hctime_leg->AddEntry(data_ep_ctime_real, Form("Reals   : %.3f", reals),"f");
  hctime_leg->AddEntry(data_ep_ctime_rand, Form("Randoms : %.3f", rands),"f");

  hctime_leg->SetBorderSize(0);
  hctime_leg->SetTextSize(0.05);
  hctime_leg->Draw();
  
  c1->Print(Form("deut_output_%s_%d.pdf", ana_type.Data(), run));
  c1->Clear();
  }
  
  // ------ INVARIANT MASS ------

  c1->cd();
  gPad->SetLogy();
  nbins = data_W_total->GetNbinsX();  //Get total number of bins (excluding overflow) (same for total, reals randoms of same histo)

  data_W_total->GetYaxis()->SetRangeUser(0.5, data_W_total->GetMaximum()+1.e5);
  data_W_total->Draw("histE0");   
  if(ana_cut!="heep_singles"){
    data_W_real->Draw("sameshistE0");   
    data_W_rand->Draw("sameshistE0");   
  }

  total = data_W_total->IntegralAndError(1, nbins, total_err);
  reals = data_W_real->IntegralAndError(1, nbins, reals_err);
  rands = data_W_rand->IntegralAndError(1, nbins, rands_err);
  
  hW_leg->AddEntry(data_W_total, Form("Total   : %.3f", total),"f");
  if(ana_cut!="heep_singles"){
    hW_leg->AddEntry(data_W_real,  Form("Reals   : %.3f", reals),"f");
    hW_leg->AddEntry(data_W_rand,  Form("Randoms : %.3f", rands),"f");
  }
  hW_leg->SetBorderSize(0);
  hW_leg->SetTextSize(0.05);
  hW_leg->Draw();

  c1->Print(Form("deut_output_%s_%d.pdf", ana_type.Data(), run));
  c1->Clear();
  
  // ------ MISSING MASS ------
  if((ana_cut!="heep_singles") || (ana_cut!="MF") || (ana_cut!="SRC")) {
  c1->cd();
  gPad->SetLogy();
  nbins = data_MM_total->GetNbinsX();  //Get total number of bins (excluding overflow) (same for total, reals randoms of same histo)

  data_MM_total->GetYaxis()->SetRangeUser(0.5, data_MM_total->GetMaximum()+1.e5);
  data_MM_total->Draw("histE0");   
  data_MM_real->Draw("sameshistE0");   
  data_MM_rand->Draw("sameshistE0");   
  
  total = data_MM_total->IntegralAndError(1, nbins, total_err);
  reals = data_MM_real->IntegralAndError(1, nbins, reals_err);
  rands = data_MM_rand->IntegralAndError(1, nbins, rands_err);
  
  hMM_leg->AddEntry(data_MM_total, Form("Total   : %.3f", total),"f");
  hMM_leg->AddEntry(data_MM_real,  Form("Reals   : %.3f", reals),"f");
  hMM_leg->AddEntry(data_MM_rand,  Form("Randoms : %.3f", rands),"f");

  hMM_leg->SetBorderSize(0);
  hMM_leg->SetTextSize(0.05);
  hMM_leg->Draw();
      
  c1->Print(Form("deut_output_%s_%d.pdf", ana_type.Data(), run));
  c1->Clear();
    
  // ------ MISSING MOMENTUM ------
  c1->cd();
  gPad->SetLogy();
  nbins = data_Pm_total->GetNbinsX();  //Get total number of bins (excluding overflow) (same for total, reals randoms of same histo)

  data_Pm_total->GetYaxis()->SetRangeUser(0.5, data_Pm_total->GetMaximum()+1.e5);
  data_Pm_total->Draw("histE0");   
  data_Pm_real->Draw("sameshistE0");   
  data_Pm_rand->Draw("sameshistE0");   
  
  total = data_Pm_total->IntegralAndError(1, nbins, total_err);
  reals = data_Pm_real->IntegralAndError(1, nbins, reals_err);
  rands = data_Pm_rand->IntegralAndError(1, nbins, rands_err);
  
  hPm_leg->AddEntry(data_Pm_total, Form("Total   : %.3f", total),"f");
  hPm_leg->AddEntry(data_Pm_real,  Form("Reals   : %.3f", reals),"f");
  hPm_leg->AddEntry(data_Pm_rand,  Form("Randoms : %.3f", rands),"f");

  hPm_leg->SetBorderSize(0);  
  hPm_leg->SetTextSize(0.05);
  hPm_leg->Draw();

  c1->Print(Form("deut_output_%s_%d.pdf", ana_type.Data(), run));
  c1->Clear();
  
  // ------ MISSING ENERGY ------
  c1->cd();
  gPad->SetLogy();
  nbins = data_Em_total->GetNbinsX();  //Get total number of bins (excluding overflow) (same for total, reals randoms of same histo)

  data_Em_total->GetYaxis()->SetRangeUser(0.5, data_Em_total->GetMaximum()+1.e5);
  data_Em_total->Draw("histE0");   
  data_Em_real->Draw("sameshistE0");   
  data_Em_rand->Draw("sameshistE0");   
  
  total = data_Em_total->IntegralAndError(1, nbins, total_err);
  reals = data_Em_real->IntegralAndError(1, nbins, reals_err);
  rands = data_Em_rand->IntegralAndError(1, nbins, rands_err);
  
  hEm_leg->AddEntry(data_Em_total, Form("Total   : %.3f", total),"f");
  hEm_leg->AddEntry(data_Em_real,  Form("Reals   : %.3f", reals),"f");
  hEm_leg->AddEntry(data_Em_rand,  Form("Randoms : %.3f", rands),"f");

  hEm_leg->SetBorderSize(0);
  hEm_leg->SetTextSize(0.05);
  hEm_leg->Draw();
  
  c1->Print(Form("deut_output_%s_%d.pdf", ana_type.Data(), run));
  c1->Clear();
  }

  if(tgt_type!="LH2"){
    // ------ 2D Nuclear Missing Energy vs. Pm ------------------------------------------
    // NOTE: Em_nuc = nu - Tp -T_{A-1}  Em_src =  nu - Tp - T_{nucleon}
    //       Em_src, is missing energy assuming kinetic energy of spectator SRC nucleon (it is used to determine a cut on to clean bkg)
    
    c1->Divide(1,2);
    
    c1->cd(1);
    //gPad->SetLogz();
    data_Em_nuc_vs_Pm->Draw("colz");
    
    c1->cd(2);
    //gPad->SetLogz();
    data_Em_src_vs_Pm->Draw("colz");
    
    c1->Print(Form("deut_output_%s_%d.pdf", ana_type.Data(), run));
    c1->Clear();
   
  //-----------------------------------------------------------------------------------
  }

  //-----------------PLOT KINEMATICS SIMC/Data comparison---------------

   //Set Legend
   auto leg_Q2 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_nu = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_W = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_thq = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_xbj = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_the = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_kf = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_Em = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_MM = new TLegend(0.1,0.8,0.28,0.9);


   c1->Divide(3,3);
   
   c1->cd(1);
   data_Q2->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
   data_Q2->GetXaxis()->CenterTitle();

   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_Q2->DrawNormalized();
     data_Q2->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_Q2->Draw();
     data_Q2->Draw("sameshistE0");
   }

   leg_Q2->AddEntry(data_Q2,"Data", "f");
   if(simc_exist) leg_Q2->AddEntry(simc_Q2,"SIMC");
   leg_Q2->Draw();
     
   c1->cd(2);
   data_nu->GetXaxis()->SetTitle("Energy Transfer, #nu [GeV]");
   data_nu->GetXaxis()->CenterTitle();  
 
   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_nu->DrawNormalized();
     data_nu->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_nu->Draw();
     data_nu->Draw("sameshistE0");
   }
   leg_nu->AddEntry(data_nu,"Data", "f");
   if(simc_exist) leg_nu->AddEntry(simc_nu,"SIMC");
   leg_nu->Draw();

   c1->cd(3);
   data_W->GetXaxis()->SetTitle("Invariant Mass, W [GeV]");
   data_W->GetXaxis()->CenterTitle();

   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_W->DrawNormalized();
     data_W->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_W->Draw();
     data_W->Draw("sameshistE0");
   }
   data_W->GetYaxis()->SetRangeUser(0., data_W->GetBinContent(data_W->GetMaximumBin())+300.);
   leg_W->AddEntry(data_W,"Data", "f");
   if(simc_exist) leg_W->AddEntry(simc_W,"SIMC");
   leg_W->Draw();

   c1->cd(4);
   data_thq->GetXaxis()->SetTitle("q-vector Angle, #theta_{q} [deg]");
   data_thq->GetXaxis()->CenterTitle();
 
   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_thq->DrawNormalized();
     data_thq->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_thq->Draw();
     data_thq->Draw("sameshistE0");
   }
   leg_thq->AddEntry(data_thq,"Data", "f");
   if(simc_exist) leg_thq->AddEntry(simc_thq,"SIMC");
   leg_thq->Draw();

   c1->cd(5);

   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_xbj->DrawNormalized();
     data_xbj->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_xbj->Draw();
     data_xbj->Draw("sameshistE0");
   }
   leg_xbj->AddEntry(data_xbj,"Data","f");
   if(simc_exist) leg_xbj->AddEntry(simc_xbj,"SIMC");
   leg_xbj->Draw();

   c1->cd(6);
   data_the->GetXaxis()->SetTitle("Electron Scatt. Angle, #theta_{e} [deg]");
   data_the->GetXaxis()->CenterTitle();

   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_the->DrawNormalized();
     data_the->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_the->Draw();
     data_the->Draw("sameshistE0");
   }
   leg_the->AddEntry(data_the,"Data","f");
   if(simc_exist) leg_the->AddEntry(simc_the,"SIMC");
   leg_the->Draw();

   c1->cd(7);
   data_kf->GetXaxis()->SetTitle("Electron Final Momentum, k_{f} [GeV/c] ");
   data_kf->GetXaxis()->CenterTitle();   
   if(simc_exist) simc_kf->Draw();
   data_kf->Draw("sameshistE0");
   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_kf->DrawNormalized();
     data_kf->DrawNormalized("sameshistE0");
   }
   leg_kf->AddEntry(data_kf,"Data","f");
   if(simc_exist) leg_kf->AddEntry(simc_kf,"SIMC");
   leg_kf->Draw();

   c1->cd(8);
   data_Em->GetXaxis()->SetTitle("Missing Energy, E_{m} [GeV/c] ");
   data_Em->GetXaxis()->CenterTitle();   

   //Draw Normalized?
   if(draw_norm){
     if((tgt_type=!"LD2") && (ana_cut=="SRC")){
       if(simc_exist) simc_Em->DrawNormalized("hist");
     }
     data_Em->DrawNormalized("sameshistE0");
     
   }
   else{
     if((tgt_type=!"LD2") && (ana_cut=="SRC")){
       if(simc_exist) simc_Em->Draw("hist");
     }
     data_Em->Draw("sameshistE0");
   }
   data_Em->GetYaxis()->SetRangeUser(0., data_Em->GetBinContent(data_Em->GetMaximumBin())+300.);
   leg_Em->AddEntry(data_Em,"Data","f");
   if((tgt_type=!"LD2") && (ana_cut=="SRC")){
     if(simc_exist) leg_Em->AddEntry(simc_Em,"SIMC");
   }
   leg_Em->Draw();

   c1->cd(9);
   data_MM->GetXaxis()->SetTitle("Missing Mass, MM [GeV] ");
   data_MM->GetXaxis()->CenterTitle();   
  
   //Draw Normalized?
   if(draw_norm){
     if((tgt_type=!"LD2") && (ana_cut=="SRC")){
       if(simc_exist) simc_MM->DrawNormalized("hist");
     }
     data_MM->DrawNormalized("sameshistE0");
   }
   else{
     if((tgt_type=!"LD2") && (ana_cut=="SRC")){
       if(simc_exist) simc_MM->Draw("hist");
     }
     data_MM->Draw("sameshistE0");
   }
   data_MM->GetYaxis()->SetRangeUser(0., data_MM->GetBinContent(data_MM->GetMaximumBin())+300.);
   leg_MM->AddEntry(data_MM,"Data","f");
   if((tgt_type=!"LD2") && (ana_cut=="SRC")){
     if(simc_exist) leg_MM->AddEntry(simc_MM,"SIMC");
   }
   leg_MM->Draw();


   c1->Print(Form("deut_output_%s_%d.pdf", ana_type.Data(), run));
   c1->Clear();                                                              

   //---------- PLOT ADDITIONAL KINEMATICS----------
   
   auto leg_Pm = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_Pf = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_thp = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_q = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_thxq = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_thrq = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_Pmx = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_Pmy = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_Pmz = new TLegend(0.1,0.8,0.28,0.9);

   c1->Divide(3,3);
   c1->cd(1);
   data_Pm->GetXaxis()->SetTitle("Missing Momentum, P_{miss} [GeV]");
   data_Pm->GetXaxis()->CenterTitle();

   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_Pm->DrawNormalized();
     data_Pm->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_Pm->Draw();
     data_Pm->Draw("sameshistE0");
   }
   data_Pm->GetYaxis()->SetRangeUser(0., data_Pm->GetBinContent(data_Pm->GetMaximumBin())+300.);
   leg_Pm->AddEntry(data_Pm,"Data", "f");
   if(simc_exist) leg_Pm->AddEntry(simc_Pm,"SIMC");
   leg_Pm->Draw();

   c1->cd(2);
   data_Pf->GetXaxis()->SetTitle("Proton Momentum, P_{p} [GeV]");
   data_Pf->GetXaxis()->CenterTitle();
 
   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_Pf->DrawNormalized();
     data_Pf->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_Pf->Draw();
     data_Pf->Draw("sameshistE0");
   }
   leg_Pf->AddEntry(data_Pf,"Data", "f");
   if(simc_exist) leg_Pf->AddEntry(simc_Pf,"SIMC");
   leg_Pf->Draw();

   c1->cd(3);
   data_thx->GetXaxis()->SetTitle("Proton Scatt. Angle, #theta_{p} [deg]");
   data_thx->GetXaxis()->CenterTitle();
 
   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_thx->DrawNormalized();
     data_thx->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_thx->Draw();
     data_thx->Draw("sameshistE0");
   }
   leg_thp->AddEntry(data_thx,"Data", "f");
   if(simc_exist) leg_thp->AddEntry(simc_thx,"SIMC");
   leg_thp->Draw();

   c1->cd(4);
   data_q->GetXaxis()->SetTitle("q-Vector Magnitude, |q| [GeV]");
   data_q->GetXaxis()->CenterTitle();
  
   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_q->DrawNormalized();
     data_q->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_q->Draw();
     data_q->Draw("sameshistE0");
   }
   leg_q->AddEntry(data_q,"Data", "f");
   if(simc_exist) leg_q->AddEntry(simc_q,"SIMC");
   leg_q->Draw();


   c1->cd(5);
   data_thxq->GetXaxis()->SetTitle("Proton-qVec. Angle, #theta_{pq} [deg]");
   data_thxq->GetXaxis()->CenterTitle();
 
   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_thxq->DrawNormalized();
     data_thxq->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_thxq->Draw();
     data_thxq->Draw("sameshistE0");
   }
   leg_thxq->AddEntry(data_thxq,"Data", "f");
   if(simc_exist) leg_thxq->AddEntry(simc_thxq,"SIMC");
   leg_thxq->Draw();

   c1->cd(6);
   data_thrq->GetXaxis()->SetTitle("Recoil-qVec. Angle, #theta_{rq} [deg]");
   data_thrq->GetXaxis()->CenterTitle();

   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_thrq->DrawNormalized();
     data_thrq->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_thrq->Draw();
     data_thrq->Draw("sameshistE0");
   }
   leg_thrq->AddEntry(data_thrq,"Data", "f");
   if(simc_exist) leg_thrq->AddEntry(simc_thrq,"SIMC");
   leg_thrq->Draw();
   
   
   c1->cd(7);
   data_Pmx->GetXaxis()->SetTitle("Missing Momentum X-comp., Pm_{x} [GeV]");
   data_Pmx->GetXaxis()->CenterTitle();
  
   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_Pmx->DrawNormalized();
     data_Pmx->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_Pmx->Draw();
     data_Pmx->Draw("sameshistE0");
   }
   leg_Pmx->AddEntry(data_Pmx,"Data", "f");
   if(simc_exist) leg_Pmx->AddEntry(simc_Pmx,"SIMC");
   leg_Pmx->Draw();

   c1->cd(8);
   data_Pmy->GetXaxis()->SetTitle("Missing Momentum Y-comp., Pm_{y} [GeV]");
   data_Pmy->GetXaxis()->CenterTitle();

   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_Pmy->DrawNormalized();
     data_Pmy->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_Pmy->Draw();
     data_Pmy->Draw("sameshistE0");
   }
   leg_Pmy->AddEntry(data_Pmy,"Data", "f");
   if(simc_exist) leg_Pmy->AddEntry(simc_Pmy,"SIMC");
   leg_Pmy->Draw();

   c1->cd(9);
   data_Pmz->GetXaxis()->SetTitle("Missing Momentum Z-comp., Pm_{z} [GeV]");
   data_Pmz->GetXaxis()->CenterTitle();

   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_Pmz->DrawNormalized();
     data_Pmz->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_Pmz->Draw();
     data_Pmz->Draw("sameshistE0");
   }
   leg_Pmz->AddEntry(data_Pmz,"Data", "f");
   if(simc_exist) leg_Pmz->AddEntry(simc_Pmz,"SIMC");
   leg_Pmz->Draw();
                                                                    

   c1->Print(Form("deut_output_%s_%d.pdf", ana_type.Data(), run));
   c1->Clear();
  
   
   //-----------------PLOT ELECTRON ARM FOCAL PLANE  Variables SIMC/Data comparison-----------------------

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
   
   c1->Divide(2,2);
   
   c1->cd(1);
 
   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_exfp->DrawNormalized();
     data_exfp->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_exfp->Draw();
     data_exfp->Draw("sameshistE0");
   }
   leg13->AddEntry(data_exfp,"Data","f");
   if(simc_exist) leg13->AddEntry(simc_exfp,"SIMC");
   leg13->Draw();
   
   c1->cd(2);
 
   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_eyfp->DrawNormalized();
     data_eyfp->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_eyfp->Draw();
     data_eyfp->Draw("sameshistE0");
   }
   leg14->AddEntry(data_eyfp,"Data", "f");
   if(simc_exist) leg14->AddEntry(simc_eyfp,"SIMC");
   leg14->Draw();

   c1->cd(3);

   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_expfp->DrawNormalized();
     data_expfp->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_expfp->Draw();
     data_expfp->Draw("sameshistE0");
   }
   leg15->AddEntry(data_expfp,"Data", "f");
   if(simc_exist) leg15->AddEntry(simc_expfp,"SIMC");
   leg15->Draw();
     
   c1->cd(4);
  
   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_eypfp->DrawNormalized();
     data_eypfp->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_eypfp->Draw();
     data_eypfp->Draw("sameshistE0");
   }
   leg16->AddEntry(data_eypfp,"Data", "f");
   if(simc_exist) leg16->AddEntry(simc_eypfp,"SIMC");
   leg16->Draw();

   c1->Print(Form("deut_output_%s_%d.pdf", ana_type.Data(), run));
   c1->Clear();
                                                                             

   //----------------------------------------------------------- 

   //-----------------PLOT HADRON ARM FOCAL PLANE  Variables SIMC/Data comparison-----------------------

   //Set Legend
   auto hfp_l1 = new TLegend(0.1,0.8,0.28,0.9);
   auto hfp_l2 = new TLegend(0.1,0.8,0.28,0.9);
   auto hfp_l3 = new TLegend(0.1,0.8,0.28,0.9);
   auto hfp_l4 = new TLegend(0.1,0.8,0.28,0.9);

   c1->Divide(2,2);

   c1->cd(1);
  
   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_hxfp->DrawNormalized();
     data_hxfp->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_hxfp->Draw();
     data_hxfp->Draw("sameshistE0");
   }
   hfp_l1->AddEntry(data_hxfp,"Data","f");
   if(simc_exist) hfp_l1->AddEntry(simc_hxfp,"SIMC");
   hfp_l1->Draw();
   
   c1->cd(2);
  
   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_hyfp->DrawNormalized();
     data_hyfp->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_hyfp->Draw();
     data_hyfp->Draw("sameshistE0");
   }
   hfp_l2->AddEntry(data_hyfp,"Data", "f");
   if(simc_exist) hfp_l2->AddEntry(simc_hyfp,"SIMC");
   hfp_l2->Draw();

   c1->cd(3);
 
   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_hxpfp->DrawNormalized();
     data_hxpfp->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_hxpfp->Draw();
     data_hxpfp->Draw("sameshistE0");
   }
   hfp_l3->AddEntry(data_hxpfp,"Data", "f");
   if(simc_exist) hfp_l3->AddEntry(simc_hxpfp,"SIMC");
   hfp_l3->Draw();
     
   c1->cd(4);
 
   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_hypfp->DrawNormalized();
     data_hypfp->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_hypfp->Draw();
     data_hypfp->Draw("sameshistE0");
   }
   hfp_l4->AddEntry(data_hypfp,"Data", "f");
   if(simc_exist) hfp_l4->AddEntry(simc_hypfp,"SIMC");
   hfp_l4->Draw();

   c1->Print(Form("deut_output_%s_%d.pdf", ana_type.Data(), run));
   c1->Clear();
                                                                                  

   //----------------------------------------------------------- 

   
   //-----------------PLOT ELECTRON Target Reconstructed Variables SIMC/Data comparison-----------------------

   //Set Legend
   auto leg5 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg6 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg7 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg8 = new TLegend(0.1,0.8,0.28,0.9);

   c1->Divide(2,2);

   c1->cd(1);
 
   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_eytar->DrawNormalized();
     data_eytar->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_eytar->Draw();
     data_eytar->Draw("sameshistE0");
   }
   leg5->AddEntry(data_eytar,"Data","f");
   if(simc_exist) leg5->AddEntry(simc_eytar,"SIMC");
   leg5->Draw();

   c1->cd(2);
 
   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_exptar->DrawNormalized();
     data_exptar->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_exptar->Draw();
     data_exptar->Draw("sameshistE0");
   }
   leg5->AddEntry(data_exptar,"Data", "f");
   if(simc_exist) leg5->AddEntry(simc_exptar,"SIMC");
   leg5->Draw();

   c1->cd(3);
 
   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_eyptar->DrawNormalized();
     data_eyptar->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_eyptar->Draw();
     data_eyptar->Draw("sameshistE0");
   }
   leg7->AddEntry(data_eyptar,"Data", "f");
   if(simc_exist) leg7->AddEntry(simc_eyptar,"SIMC");
   leg7->Draw();
     
   c1->cd(4);
 
   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_edelta->DrawNormalized();
     data_edelta->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_edelta->Draw();
     data_edelta->Draw("sameshistE0");
   }
   leg8->AddEntry(data_edelta,"Data", "f");
   if(simc_exist) leg8->AddEntry(simc_edelta,"SIMC");
   leg8->Draw();
   
   c1->Print(Form("deut_output_%s_%d.pdf", ana_type.Data(), run));
   c1->Clear();


   //------------------------------------------------------------------------------

   //-----------------PLOT Target Reconstructed Variables SIMC/Data comparison-----------------------
 
   //Set Legend
   auto htr_l1 = new TLegend(0.1,0.8,0.28,0.9);
   auto htr_l2 = new TLegend(0.1,0.8,0.28,0.9);
   auto htr_l3 = new TLegend(0.1,0.8,0.28,0.9);
   auto htr_l4 = new TLegend(0.1,0.8,0.28,0.9);

   c1->Divide(2,2);

   c1->cd(1);
  
   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_hytar->DrawNormalized();
     data_hytar->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_hytar->Draw();
     data_hytar->Draw("sameshistE0");
   }
   htr_l1->AddEntry(data_hytar,"Data","f");
   if(simc_exist) htr_l1->AddEntry(simc_hytar,"SIMC");
   htr_l1->Draw();

   c1->cd(2);
   
   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_hxptar->DrawNormalized();
     data_hxptar->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_hxptar->Draw();
     data_hxptar->Draw("sameshistE0");
   }
   htr_l2->AddEntry(data_hxptar,"Data", "f");
   if(simc_exist) htr_l2->AddEntry(simc_hxptar,"SIMC");
   htr_l2->Draw();

   c1->cd(3);
  
   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_hyptar->DrawNormalized();
     data_hyptar->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_hyptar->Draw();
     data_hyptar->Draw("sameshistE0");
   }
   htr_l3->AddEntry(data_hyptar,"Data", "f");
   if(simc_exist) htr_l3->AddEntry(simc_hyptar,"SIMC");
   htr_l3->Draw();
     
   c1->cd(4);
 
   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_hdelta->DrawNormalized();
     data_hdelta->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_hdelta->Draw();
     data_hdelta->Draw("sameshistE0");
   }
   htr_l4->AddEntry(data_hdelta,"Data", "f");
   if(simc_exist) htr_l4->AddEntry(simc_hdelta,"SIMC");
   htr_l4->Draw();

   c1->Print(Form("deut_output_%s_%d.pdf", ana_type.Data(), run));
   c1->Clear();
   

   //------------------------------------------------------------------------------

 //-----------------PLOT TARGET VERTEX Variables SIMC/Data comparison-----------------------

  //Set Legend
   auto leghxt = new TLegend(0.1,0.8,0.28,0.9);                          
   auto leghyt = new TLegend(0.1,0.8,0.28,0.9);  
   auto leghzt = new TLegend(0.1,0.8,0.28,0.9);                                                          
                                                                                                                                                          
   auto legpxt = new TLegend(0.1,0.8,0.28,0.9);                              
   auto legpyt = new TLegend(0.1,0.8,0.28,0.9);                                                                       
   auto legpzt = new TLegend(0.1,0.8,0.28,0.9);   

   c1->Divide(3,2);
   
   c1->cd(1);
 
   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_xtarH->DrawNormalized("hist");
     data_xtarH->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_xtarP->Draw("hist");
     data_xtarH->Draw("sameshistE0");
   }
   leghxt->AddEntry(data_xtarH,"Data","f");
   if(simc_exist) leghxt->AddEntry(simc_xtarH,"SIMC");
   leghxt->Draw();
  
   c1->cd(2);
  
   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_ytarH->DrawNormalized("hist");
     data_ytarH->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_ytarH->Draw("hist");
     data_ytarH->Draw("sameshistE0");
   }
   leghyt->AddEntry(data_ytarH,"Data","f");
   if(simc_exist) leghyt->AddEntry(simc_ytarH,"SIMC");
   leghyt->Draw();

   c1->cd(3);
 
   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_ztarH->DrawNormalized("hist");
     data_ztarH->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_ztarH->Draw("hist");
     data_ztarH->Draw("sameshistE0");
   }
   leghzt->AddEntry(data_ztarH,"Data","f");
   if(simc_exist) leghzt->AddEntry(simc_ztarH,"SIMC");
   leghzt->Draw();   

   c1->cd(4);

   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_xtarP->DrawNormalized("hist");
     data_xtarP->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_xtarP->Draw("hist");
     data_xtarP->Draw("sameshistE0");
   }
   legpxt->AddEntry(data_xtarP,"Data","f");
   if(simc_exist) legpxt->AddEntry(simc_xtarP,"SIMC");
   legpxt->Draw();
  
   c1->cd(5);

   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_ytarP->DrawNormalized("hist");
     data_ytarP->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_ytarP->Draw("hist");
     data_ytarP->Draw("sameshistE0");
   }
   legpyt->AddEntry(data_ytarP,"Data","f");
   if(simc_exist) legpyt->AddEntry(simc_ytarP,"SIMC");
   legpyt->Draw();

   c1->cd(6);
 
   //Draw Normalized?
   if(draw_norm){
     if(simc_exist) simc_ztarP->DrawNormalized("hist");
     data_ztarP->DrawNormalized("sameshistE0");
   }
   else{
     if(simc_exist) simc_ztarP->Draw("hist");
     data_ztarP->Draw("sameshistE0");
   }
   legpzt->AddEntry(data_ztarP,"Data","f");
   if(simc_exist) legpzt->AddEntry(simc_ztarP,"SIMC");
   legpzt->Draw();

   c1->Print(Form("deut_output_%s_%d.pdf", ana_type.Data(), run));
   c1->Clear();

   //--------------------------------------------------------------------
   
  

   // Complete writing out multi-page .pdf
   c1->Print(Form("deut_output_%s_%d.pdf]", ana_type.Data(), run));
   
   gSystem->Exec(Form("mv deut_output_%s_%d.pdf %s", ana_type.Data(), run, outPDF.Data()));
   gSystem->Exec(Form("evince %s", outPDF.Data()));
   //gSystem->Exec(Form("open deut_output_%s_%d.pdf", ana_type.Data(), run));
      
}
