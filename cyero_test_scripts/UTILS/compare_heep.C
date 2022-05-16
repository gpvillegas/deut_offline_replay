//Script to make comparison between SIMC and Commissioning Data from HallC Spring 2018
//Compare Target Reconstruction/FOCAL PLANE/ Kinematics Variables

void compare_heep(int run)
{

  gROOT->SetBatch(kTRUE);  
  gStyle->SetOptStat(1001111);
  //TString simc_filename =  "weighted_ep_coin_simc_1854.root"; //"ep_coin_simc_1929.root";
  
  //Pre-defined SIMC/data root file names containing histogram object to comapare
  TString simc_filename =  Form("../heep_simc_histos_%d_rad.root", run);                      
  TString data_filename = Form("../heep_data_histos_%d_combined.root",run); 

  //Where to store plots
  string plots_dir = "./";
  string plots_path;
  
  //Open SIMC/data ROOT files;
  TFile *simc_file = new TFile(simc_filename);
  TFile *data_file = new TFile(data_filename);

  //---------------Target ----------------
  //Define SIMC histos ('h'-->hadron arm,  'e'-->electron arm)
  
  TH1F *simc_xtar =  0;
  TH1F *simc_ytarH =  0;
  TH1F *simc_ztarH =  0;

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
  TH1F *simc_omega =  0;
  TH1F *simc_W =  0;
  TH1F *simc_thq = 0;

  TH1F *simc_xbj = 0;
  TH1F *simc_th_elec = 0;                                  
  TH1F *simc_kf = 0;  
  TH1F *simc_emiss = 0;

  //Kinematics 2
  TH1F *simc_Pm = 0;
  TH1F *simc_Pf = 0;
  TH1F *simc_th_prot = 0;
  TH1F *simc_q = 0;    //q-vector magnitude
  TH1F *simc_thpq = 0;
  TH1F *simc_Pmx = 0;
  TH1F *simc_Pmy = 0;
  TH1F *simc_Pmz = 0;

  //Define data histos
  TH1F *data_Q2 =  0;
  TH1F *data_omega =  0;
  TH1F *data_W =  0;
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
  TH1F *data_Pmx = 0;
  TH1F *data_Pmy = 0;
  TH1F *data_Pmz = 0;

  //---------------------------------------------------------------

 //change to simc_file
  simc_file->cd();

  
  //----------Get Target Histograms------------------
  //Get Histogram objects from SIMC rootfile
  simc_file->GetObject("H_hx_tar", simc_xtar);

  simc_file->GetObject("H_hy_tar", simc_ytarH);
  simc_file->GetObject("H_hz_tar", simc_ztarH);

  simc_file->GetObject("H_ey_tar", simc_ytarP);  
  simc_file->GetObject("H_ez_tar", simc_ztarP); 

  //Set SIMC Histo Aesthetics
  simc_xtar->SetLineColor(kRed);
  simc_xtar->SetLineWidth(2);
 
  simc_ytarH->SetLineColor(kRed);
  simc_ytarH->SetLineWidth(2);
  simc_ztarH->SetLineColor(kRed);
  simc_ztarH->SetLineWidth(2);

  simc_ytarP->SetLineColor(kRed);          
  simc_ytarP->SetLineWidth(2);           
  simc_ztarP->SetLineColor(kRed);                        
  simc_ztarP->SetLineWidth(2);  
  
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

 //change to simc_file
  simc_file->cd();

  //----------Get Target Reconstructed Histograms------------------
  //Get Histogram objects from SIMC rootfile
  simc_file->GetObject("H_eytar", simc_eytar);
  simc_file->GetObject("H_exptar", simc_exptar);
  simc_file->GetObject("H_eyptar", simc_eyptar);
  simc_file->GetObject("H_edelta", simc_edelta);

  simc_file->GetObject("H_hytar", simc_hytar);
  simc_file->GetObject("H_hxptar", simc_hxptar);
  simc_file->GetObject("H_hyptar", simc_hyptar);
  simc_file->GetObject("H_hdelta", simc_hdelta);

  simc_file->GetObject("H_Pmx_Lab", simc_Pmx);
  simc_file->GetObject("H_Pmy_Lab", simc_Pmy);
  simc_file->GetObject("H_Pmz_Lab", simc_Pmz);

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

  simc_Pmx->SetLineColor(kRed);
  simc_Pmx->SetLineWidth(2);
  simc_Pmy->SetLineColor(kRed);
  simc_Pmy->SetLineWidth(2);
  simc_Pmz->SetLineColor(kRed);
  simc_Pmz->SetLineWidth(2);

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

   //change to simc_file
  simc_file->cd();

  //Get Histogram objects from SIMC rootfile
  simc_file->GetObject("H_exfp", simc_exfp);
  simc_file->GetObject("H_eyfp", simc_eyfp);
  simc_file->GetObject("H_expfp", simc_expfp);
  simc_file->GetObject("H_eypfp", simc_eypfp);

  simc_file->GetObject("H_hxfp", simc_hxfp);
  simc_file->GetObject("H_hyfp", simc_hyfp);
  simc_file->GetObject("H_hxpfp", simc_hxpfp);
  simc_file->GetObject("H_hypfp", simc_hypfp);
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

   //change to simc_file
  simc_file->cd();

  //Get Histogram objects from SIMC rootfile
  simc_file->GetObject("H_Q2", simc_Q2);
  simc_file->GetObject("H_omega", simc_omega);
  simc_file->GetObject("H_W", simc_W);
  simc_file->GetObject("H_theta_q", simc_thq);

  simc_file->GetObject("H_xbj", simc_xbj);
  simc_file->GetObject("H_theta_elec", simc_th_elec);
  simc_file->GetObject("H_kf", simc_kf);
  simc_file->GetObject("H_Emiss", simc_emiss);

  simc_file->GetObject("H_Pm", simc_Pm);
  simc_file->GetObject("H_Pf", simc_Pf);
  simc_file->GetObject("H_theta_prot", simc_th_prot);
  simc_file->GetObject("H_q", simc_q);
  simc_file->GetObject("H_theta_pq", simc_thpq);

  
  
  //Set SIMC Histo Aesthetics
  simc_Q2->SetLineColor(kRed);
  simc_Q2->SetLineWidth(2);
  simc_omega->SetLineColor(kRed);
  simc_omega->SetLineWidth(2);
  simc_W->SetLineColor(kRed);
  simc_W->SetLineWidth(2);
  simc_thq->SetLineColor(kRed);
  simc_thq->SetLineWidth(2);
  
  simc_xbj->SetLineColor(kRed);
  simc_xbj->SetLineWidth(2);
  simc_th_elec->SetLineColor(kRed);
  simc_th_elec->SetLineWidth(2);
  simc_kf->SetLineColor(kRed);
  simc_kf->SetLineWidth(2);
  simc_emiss->SetLineColor(kRed);
  simc_emiss->SetLineWidth(2);
  
  simc_Pm->SetLineColor(kRed);
  simc_Pm->SetLineWidth(2);
  simc_Pf->SetLineColor(kRed);
  simc_Pf->SetLineWidth(2);
  simc_th_prot->SetLineColor(kRed);
  simc_th_prot->SetLineWidth(2);
  simc_q->SetLineColor(kRed);
  simc_q->SetLineWidth(2);
  simc_thpq->SetLineColor(kRed);
  simc_thpq->SetLineWidth(2);

  //change to data_file
  data_file->cd();
  
  //Get Histogram objects from data rootfile
  data_file->GetObject("H_Q2", data_Q2);
  data_file->GetObject("H_omega", data_omega);
  data_file->GetObject("H_W", data_W);
  data_file->GetObject("H_theta_q", data_thq);

  
  data_file->GetObject("H_xbj", data_xbj);
  data_file->GetObject("H_theta_elec", data_th_elec);
  data_file->GetObject("H_kf", data_kf);
  data_file->GetObject("H_Emiss", data_emiss);

  data_file->GetObject("H_Pm", data_Pm);
  data_file->GetObject("H_Pf", data_Pf);
  data_file->GetObject("H_theta_prot", data_th_prot);
  data_file->GetObject("H_q", data_q);
  data_file->GetObject("H_theta_pq", data_thpq);

  data_file->GetObject("H_Pmx_Lab", data_Pmx);
  data_file->GetObject("H_Pmy_Lab", data_Pmy);
  data_file->GetObject("H_Pmz_Lab", data_Pmz);

  //Set data Histo Aesthetics
  data_Q2->SetFillColorAlpha(kBlue, 0.35);
  data_Q2->SetFillStyle(3004);
  data_omega->SetFillColorAlpha(kBlue, 0.35);
  data_omega->SetFillStyle(3004);
  data_W->SetFillColorAlpha(kBlue, 0.35);
  data_W->SetFillStyle(3004);
  data_thq->SetFillColorAlpha(kBlue, 0.35);
  data_thq->SetFillStyle(3004);

  data_xbj->SetFillColorAlpha(kBlue, 0.35);
  data_xbj->SetFillStyle(3004);
  data_th_elec->SetFillColorAlpha(kBlue, 0.35);
  data_th_elec->SetFillStyle(3004);
  data_kf->SetFillColorAlpha(kBlue, 0.35);
  data_kf->SetFillStyle(3004);
  data_emiss->SetFillColorAlpha(kBlue, 0.35);
  data_emiss->SetFillStyle(3004);

  data_Pm->SetFillColorAlpha(kBlue,0.35);
  data_Pm->SetFillStyle(3004);
  data_Pf->SetFillColorAlpha(kBlue,0.35);
  data_Pf->SetFillStyle(3004);
  data_th_prot->SetFillColorAlpha(kBlue,0.35);
  data_th_prot->SetFillStyle(3004);
  data_q->SetFillColorAlpha(kBlue,0.35);
  data_q->SetFillStyle(3004);
  data_thpq->SetFillColorAlpha(kBlue,0.35);
  data_thpq->SetFillStyle(3004);

  data_Pmx->SetFillColorAlpha(kBlue,0.35);
  data_Pmx->SetFillStyle(3004);
  data_Pmy->SetFillColorAlpha(kBlue,0.35);
  data_Pmy->SetFillStyle(3004);
  data_Pmz->SetFillColorAlpha(kBlue,0.35);
  data_Pmz->SetFillStyle(3004);

  //Overlay SIMC/data plots (*** VERY IMPORTANT ***: Range and #bins must be same)


   //Set Legend
   auto leg5 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg6 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg7 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg8 = new TLegend(0.1,0.8,0.28,0.9);

  
   //-----------------PLOT Target Reconstructed Variables SIMC/Data comparison-----------------------

   //Create A Canvas to store Target Recon. variable comparisons in HADRON ARM
   
   TCanvas *c1 = new TCanvas("c1", "Electron Arm: Target Reconstruction", 2000, 1000);
   c1->Divide(2,2);

   c1->cd(1);
   simc_eytar->Draw();
   data_eytar->Draw("sameshist");
   leg5->AddEntry(data_eytar,"Data","f");
   leg5->AddEntry(simc_eytar,"SIMC");
   leg5->Draw();

   c1->cd(2);
   simc_exptar->Draw();
   data_exptar->Draw("sameshist");
   leg5->AddEntry(data_exptar,"Data", "f");
   leg5->AddEntry(simc_exptar,"SIMC");
   leg5->Draw();

   c1->cd(3);
   simc_eyptar->Draw();
   data_eyptar->Draw("sameshist");
   leg7->AddEntry(data_eyptar,"Data", "f");
   leg7->AddEntry(simc_eyptar,"SIMC");
   leg7->Draw();
     
   c1->cd(4);
   simc_edelta->Draw();
   data_edelta->Draw("sameshist");
   leg8->AddEntry(data_edelta,"Data", "f");
   leg8->AddEntry(simc_edelta,"SIMC");
   leg8->Draw();

   plots_path = plots_dir + Form("eArm_TargRecon_%d.pdf", run);
   c1->SaveAs(plots_path.c_str());

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

   TCanvas *c2 = new TCanvas("c2", "Electron Arm: Focal Plane", 2000, 1000);
   c2->Divide(2,2);

   c2->cd(1);
   simc_exfp->Draw();
   data_exfp->Draw("sameshist");
   leg13->AddEntry(data_exfp,"Data","f");
   leg13->AddEntry(simc_exfp,"SIMC");
   leg13->Draw();
   
   c2->cd(2);
   simc_eyfp->Draw();
   data_eyfp->Draw("sameshist");
   leg14->AddEntry(data_eyfp,"Data", "f");
   leg14->AddEntry(simc_eyfp,"SIMC");
   leg14->Draw();

   c2->cd(3);
   simc_expfp->Draw();
   data_expfp->Draw("sameshist");
   leg15->AddEntry(data_expfp,"Data", "f");
   leg15->AddEntry(simc_expfp,"SIMC");
   leg15->Draw();
     
   c2->cd(4);
   simc_eypfp->Draw();
   data_eypfp->Draw("sameshist");
   leg16->AddEntry(data_eypfp,"Data", "f");
   leg16->AddEntry(simc_eypfp,"SIMC");
   leg16->Draw();

   plots_path = plots_dir + Form("eArm_FocalPlane_%d.pdf", run);
   c2->SaveAs(plots_path.c_str());                                                                                   

   //----------------------------------------------------------- 
 
   
   //-----------------PLOT KINEMATICS SIMC/Data comparison---------------

//Set Legend
   auto leg19 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg20 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg21 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg22 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg23 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg24 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg25 = new TLegend(0.1,0.8,0.28,0.9);

   TCanvas *c3 = new TCanvas("c3", "Kinematics", 2000, 1000);
   c3->Divide(4,2);
   
   c3->cd(1);
   data_Q2->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
   data_Q2->GetXaxis()->CenterTitle();
   simc_Q2->Draw();
   data_Q2->Draw("sameshist");
   leg19->AddEntry(data_Q2,"Data", "f");
   leg19->AddEntry(simc_Q2,"SIMC");
   leg19->Draw();
     
   c3->cd(2);
   data_omega->GetXaxis()->SetTitle("Energy Transfer, #omega [GeV]");
   data_omega->GetXaxis()->CenterTitle();  
   simc_omega->Draw();
   data_omega->Draw("sameshist");
   leg20->AddEntry(data_omega,"Data", "f");
   leg20->AddEntry(simc_omega,"SIMC");
   leg20->Draw();

   c3->cd(3);
   data_W->GetXaxis()->SetTitle("Invariant Mass, W [GeV]");
   data_W->GetXaxis()->CenterTitle();
   simc_W->Draw();
   data_W->Draw("sameshist");
   leg21->AddEntry(data_W,"Data", "f");
   leg21->AddEntry(simc_W,"SIMC");
   leg21->Draw();

   c3->cd(4);
   data_thq->GetXaxis()->SetTitle("q-vector Angle, #theta_{q} [deg]");
   data_thq->GetXaxis()->CenterTitle();
   simc_thq->Draw();
   data_thq->Draw("sameshist");
   leg22->AddEntry(data_thq,"Data", "f");
   leg22->AddEntry(simc_thq,"SIMC");
   leg22->Draw();

   c3->cd(5);
   simc_xbj->Draw();
   data_xbj->Draw("sameshist");
   leg23->AddEntry(data_xbj,"Data","f");
   leg23->AddEntry(simc_xbj,"SIMC");
   leg23->Draw();

   c3->cd(6);
   data_th_elec->GetXaxis()->SetTitle("Electron Scatt. Angle, #theta_{e} [deg]");
   data_th_elec->GetXaxis()->CenterTitle();
   simc_th_elec->Draw();
   data_th_elec->Draw("sameshist");
   leg24->AddEntry(data_th_elec,"Data","f");
   leg24->AddEntry(simc_th_elec,"SIMC");
   leg24->Draw();

   c3->cd(7);
   data_kf->GetXaxis()->SetTitle("Electron Final Momentum, k_{f} [GeV/c] ");
   data_kf->GetXaxis()->CenterTitle();   
   simc_kf->Draw();
   data_kf->Draw("sameshist");
   leg25->AddEntry(data_kf,"Data","f");
   leg25->AddEntry(simc_kf,"SIMC");
   leg25->Draw();

   c3->cd(8);
   data_emiss->GetXaxis()->SetTitle("Missing Energy, E_{m} [GeV/c] ");
   data_emiss->GetXaxis()->CenterTitle();   
   simc_emiss->Draw("hist");
   data_emiss->Draw("sameshist");
   leg25->AddEntry(data_emiss,"Data","f");
   leg25->AddEntry(simc_emiss,"SIMC");
   leg25->Draw();

   plots_path = plots_dir + Form("Kinematics1_%d.pdf", run);
   c3->SaveAs(plots_path.c_str());                                                                   

   //Plot Additional Kinematics
   
   auto leg_Pm = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_Pf = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_thp = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_q = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_thpq = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_Pmx = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_Pmy = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_Pmz = new TLegend(0.1,0.8,0.28,0.9);


   //Create A Canvas to store kinematic variable comparisons
   TCanvas *ck2 = new TCanvas("ck2", "Kinematics-2", 2000, 1000);
   
   ck2->Divide(4,2);
   ck2->cd(1);
   data_Pm->GetXaxis()->SetTitle("Missing Momentum, P_{miss} [GeV]");
   data_Pm->GetXaxis()->CenterTitle();
   simc_Pm->Draw();
   data_Pm->Draw("sameshist");
   leg_Pm->AddEntry(data_Pm,"Data", "f");
   leg_Pm->AddEntry(simc_Pm,"SIMC");
   leg_Pm->Draw();

   ck2->cd(2);
   data_Pf->GetXaxis()->SetTitle("Proton Momentum, P_{p} [GeV]");
   data_Pf->GetXaxis()->CenterTitle();
   simc_Pf->Draw();
   data_Pf->Draw("sameshist");
   leg_Pf->AddEntry(data_Pf,"Data", "f");
   leg_Pf->AddEntry(simc_Pf,"SIMC");
   leg_Pf->Draw();

   ck2->cd(3);
   data_th_prot->GetXaxis()->SetTitle("Proton Scatt. Angle, #theta_{p} [deg]");
   data_th_prot->GetXaxis()->CenterTitle();
   simc_th_prot->Draw();
   data_th_prot->Draw("sameshist");
   leg_thp->AddEntry(data_th_prot,"Data", "f");
   leg_thp->AddEntry(simc_th_prot,"SIMC");
   leg_thp->Draw();

   ck2->cd(4);
   data_q->GetXaxis()->SetTitle("q-Vector Magnitude, |q| [GeV]");
   data_q->GetXaxis()->CenterTitle();
   simc_q->Draw();
   data_q->Draw("sameshist");
   leg_q->AddEntry(data_q,"Data", "f");
   leg_q->AddEntry(simc_q,"SIMC");
   leg_q->Draw();


   ck2->cd(5);
   data_thpq->GetXaxis()->SetTitle("Proton-qVec. Angle, #theta_{pq} [deg]");
   data_thpq->GetXaxis()->CenterTitle();
   simc_thpq->Draw();
   data_thpq->Draw("sameshist");
   leg_thpq->AddEntry(data_thpq,"Data", "f");
   leg_thpq->AddEntry(simc_thpq,"SIMC");
   leg_thpq->Draw();
  
   ck2->cd(6);
   data_Pmx->GetXaxis()->SetTitle("Missing Momentum X-comp., Pm_{x} [GeV]");
   data_Pmx->GetXaxis()->CenterTitle();
   simc_Pmx->Draw();
   data_Pmx->Draw("sameshist");
   leg_Pmx->AddEntry(data_Pmx,"Data", "f");
   leg_Pmx->AddEntry(simc_Pmx,"SIMC");
   leg_Pmx->Draw();

   ck2->cd(7);
   data_Pmy->GetXaxis()->SetTitle("Missing Momentum Y-comp., Pm_{y} [GeV]");
   data_Pmy->GetXaxis()->CenterTitle();
   simc_Pmy->Draw();
   data_Pmy->Draw("sameshist");
   leg_Pmy->AddEntry(data_Pmy,"Data", "f");
   leg_Pmy->AddEntry(simc_Pmy,"SIMC");
   leg_Pmy->Draw();

   ck2->cd(8);
   data_Pmz->GetXaxis()->SetTitle("Missing Momentum Z-comp., Pm_{z} [GeV]");
   data_Pmz->GetXaxis()->CenterTitle();
   simc_Pmz->Draw();
   data_Pmz->Draw("sameshist");
   leg_Pmz->AddEntry(data_Pmz,"Data", "f");
   leg_Pmz->AddEntry(simc_Pmz,"SIMC");
   leg_Pmz->Draw();
   
   plots_path = plots_dir + Form("Kinematics2_%d.pdf", run);
   ck2->SaveAs(plots_path.c_str());                                                                   


 //-----------------PLOT TARGET  Variables SIMC/Data comparison-----------------------

  //Set Legend
   auto leghxt = new TLegend(0.1,0.8,0.28,0.9);                          
   auto leghyt = new TLegend(0.1,0.8,0.28,0.9);  
   auto leghzt = new TLegend(0.1,0.8,0.28,0.9);                                                          
                                                                                                                                                          
   auto legpxt = new TLegend(0.1,0.8,0.28,0.9);                              
   auto legpyt = new TLegend(0.1,0.8,0.28,0.9);                                                                       
   auto legpzt = new TLegend(0.1,0.8,0.28,0.9);   

TCanvas *c4a = new TCanvas("c4a", "HMS Target Variables", 2000, 1000);
   c4a->Divide(3,1);

   c4a->cd(1);
   simc_xtar->Draw("hist");
   data_xtarH->Draw("sameshist");
   leghxt->AddEntry(data_xtarH,"Data","f");
   leghxt->AddEntry(simc_xtar,"SIMC");
   leghxt->Draw();
  
   c4a->cd(2);
   simc_ytarH->Draw("hist");
   data_ytarH->Draw("sameshist");
   leghyt->AddEntry(data_ytarH,"Data","f");
   leghyt->AddEntry(simc_ytarH,"SIMC");
   leghyt->Draw();

   c4a->cd(3);
   simc_ztarH->Draw("hist");
   data_ztarH->Draw("sameshist");
   leghzt->AddEntry(data_ztarH,"Data","f");
   leghzt->AddEntry(simc_ztarH,"SIMC");
   leghzt->Draw();
   
   plots_path = plots_dir + Form("hArm_TargVar_%d.pdf", run);
   c4a->SaveAs(plots_path.c_str());                                                                                              

   TCanvas *c4b = new TCanvas("c4b", "SHMS Target Variables", 2000, 1000);
   c4b->Divide(3,1);

   c4b->cd(1);
   simc_xtar->Draw("hist");
   data_xtarP->Draw("sameshist");
   legpxt->AddEntry(data_xtarP,"Data","f");
   legpxt->AddEntry(simc_xtar,"SIMC");
   legpxt->Draw();
  
   c4b->cd(2);
   simc_ytarP->Draw("hist");
   data_ytarP->Draw("sameshist");
   legpyt->AddEntry(data_ytarP,"Data","f");
   legpyt->AddEntry(simc_ytarP,"SIMC");
   legpyt->Draw();

   c4b->cd(3);
   simc_ztarP->Draw("hist");
   data_ztarP->Draw("sameshist");
   legpzt->AddEntry(data_ztarP,"Data","f");
   legpzt->AddEntry(simc_ztarP,"SIMC");
   legpzt->Draw();
   
   plots_path = plots_dir + Form("pArm_TargVar_%d.pdf", run);
   c4b->SaveAs(plots_path.c_str());      
   //--------PLOT HADRON ARM QUANTITIES--------



   //--------PLOT HADRON ARM QUANTITIES--------


   
   //-----------------PLOT Target Reconstructed Variables SIMC/Data comparison-----------------------
 
   //Set Legend
   auto htr_l1 = new TLegend(0.1,0.8,0.28,0.9);
   auto htr_l2 = new TLegend(0.1,0.8,0.28,0.9);
   auto htr_l3 = new TLegend(0.1,0.8,0.28,0.9);
   auto htr_l4 = new TLegend(0.1,0.8,0.28,0.9);
   
   //Create A Canvas to store Target Recon. variable comparisons in HADRON ARM
   
   TCanvas *htr = new TCanvas("htr", "Hadron Arm: Target Reconstruction", 2000, 1000);
   htr->Divide(2,2);

   htr->cd(1);
   simc_hytar->Draw();
   data_hytar->Draw("sameshist");
   htr_l1->AddEntry(data_hytar,"Data","f");
   htr_l1->AddEntry(simc_hytar,"SIMC");
   htr_l1->Draw();

   htr->cd(2);
   simc_hxptar->Draw();
   data_hxptar->Draw("sameshist");
   htr_l2->AddEntry(data_hxptar,"Data", "f");
   htr_l2->AddEntry(simc_hxptar,"SIMC");
   htr_l2->Draw();

   htr->cd(3);
   simc_hyptar->Draw();
   data_hyptar->Draw("sameshist");
   htr_l3->AddEntry(data_hyptar,"Data", "f");
   htr_l3->AddEntry(simc_hyptar,"SIMC");
   htr_l3->Draw();
     
   htr->cd(4);
   simc_hdelta->Draw();
   data_hdelta->Draw("sameshist");
   htr_l4->AddEntry(data_hdelta,"Data", "f");
   htr_l4->AddEntry(simc_hdelta,"SIMC");
   htr_l4->Draw();
   
   plots_path = plots_dir + Form("hArm_TargRecon_%d.pdf", run);
   htr->SaveAs(plots_path.c_str());

   //------------------------------------------------------------------------------

   
   //-----------------PLOT FOCAL PLANE  Variables SIMC/Data comparison-----------------------

   //Set Legend
   auto hfp_l1 = new TLegend(0.1,0.8,0.28,0.9);
   auto hfp_l2 = new TLegend(0.1,0.8,0.28,0.9);
   auto hfp_l3 = new TLegend(0.1,0.8,0.28,0.9);
   auto hfp_l4 = new TLegend(0.1,0.8,0.28,0.9);

   TCanvas *hfp = new TCanvas("hfp", "Hadron Arm: Focal Plane", 2000, 1000);
   hfp->Divide(2,2);

   hfp->cd(1);
   simc_hxfp->Draw();
   data_hxfp->Draw("sameshist");
   hfp_l1->AddEntry(data_hxfp,"Data","f");
   hfp_l1->AddEntry(simc_hxfp,"SIMC");
   hfp_l1->Draw();
   
   hfp->cd(2);
   simc_hyfp->Draw();
   data_hyfp->Draw("sameshist");
   hfp_l2->AddEntry(data_hyfp,"Data", "f");
   hfp_l2->AddEntry(simc_hyfp,"SIMC");
   hfp_l2->Draw();

   hfp->cd(3);
   simc_hxpfp->Draw();
   data_hxpfp->Draw("sameshist");
   hfp_l3->AddEntry(data_hxpfp,"Data", "f");
   hfp_l3->AddEntry(simc_hxpfp,"SIMC");
   hfp_l3->Draw();
     
   hfp->cd(4);
   simc_hypfp->Draw();
   data_hypfp->Draw("sameshist");
   hfp_l4->AddEntry(data_hypfp,"Data", "f");
   hfp_l4->AddEntry(simc_hypfp,"SIMC");
   hfp_l4->Draw();
   
   plots_path = plots_dir + Form("hArm_FocalPlane_%d.pdf", run);
   hfp->SaveAs(plots_path.c_str());                                                                                   

   //----------------------------------------------------------- 
 
  
}
