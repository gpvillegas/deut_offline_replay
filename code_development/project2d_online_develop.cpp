//#ifndef PROJECT_2D_H
//#define PROJECT_2D_H

/*
  Author: C. Yero
  Date: Feb 15, 2023 
 */
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;



// helper function : check if a string only contain numebers or not
int is_digits(string& str)
{
    for (char ch : str) {
        int v = ch; // ASCII Val converted
        if (!(ch >= 48 && ch <= 57)) {   // The ASCII value of '0' is 48 , '9' is 57.
            
	  return 0;  // characters in string have letters
        }
	
	
    }
 
    return 1; // no charactes have letter (therefore true number)
}


// helper function : combined 2d histos for a given run range based on a .csv file,
// and looks for the histos in a pre-defined path
TH2F* combine_2dhistos(int run_min=0, int run_max=99999, TString hist2d_name="randSub_plots/H_Pm_vs_thrq_rand_sub")
{

  // cout << Form("Will skip "" runs", skip_len) << endl;
  
  //experiment online file to read run numbers (assumes .csv format, and runs to be 1st column)
  ifstream file("UTILS_DEUT/runlist/deut-2023_runlist.csv"); 
  string line;
  string token;
  int run_num;

  // define generic root file name to be read 
  string root_file_path; // generic path to where root files are written

  // define tfile to be used to read .root
  TFile *data_file = NULL;

  // define generic histogram to be combined
  TH2F *myhist2d = 0;
  TH2F *myhist2d_total = 0;

  int cnt=0; // good run counter
  
    // read line by line
  while (getline(file, line)) {

    // get token (parsed line string up to the 1st comma)
    token = line.substr(0, line.find(','));

    // check if token is a number
    if( is_digits(token) ){


      // convert token to int
      run_num = stoi(token);

      
      // define generic rootfile name to path
      root_file_path =Form("DEUT_OUTPUT/ROOT/deut_prod_LD2_deep_%d_-1_histos.root", run_num);

      // if file does not exist continue
      if(gSystem->AccessPathName(root_file_path.c_str())) continue;

      // skip runs outside range
      if( (run_num<run_min) or (run_num>run_max) ) continue;
      
      
      cout << "---> reading ROOTfile: " << root_file_path.c_str() << endl;


      // for each root file, get the desired histogram to be combined
      data_file =  new TFile(root_file_path.c_str(), "READ");

      data_file->GetObject(Form("%s", hist2d_name.Data()), myhist2d);

      cout << "data_file read properly" << endl;
      // only for 1st run, clone histogram to the total
      if(cnt==0) { myhist2d_total = (TH2F*)myhist2d->Clone(hist2d_name.Data()); }

      cout << "passed L1" << endl;
      // add subseqquent histos
      if(cnt>0) { myhist2d_total->Add(myhist2d); }
      
      // increment counter
      cnt++;
      
    } // end check if token is number
    
  } // end readlines

  cout << "passed L2" << endl;

  return myhist2d_total;

}


TH2F* get_simc_2d_histos(TString setting="pm120", TString hist_type="rad_corr_ratio"){

  // Brief: this function gets specific histogram from SIMC for
  // radiative corrections and phase space, for a given deuteron kin
  // the arguments are:  setting="pm120", "pm580", "pm800", "pm900"
  // hist_type="norad", "rad", "rad_corr_ratio" , "phase_space", 

  cout << "calling get_simc_2d_histos()" << endl;
  
  // simc no_rad / rad 2D histo
  TH2F *dummy = 0;
  TH2F *H2_Pm_vs_thrq_simc_rad = 0;
  TH2F *H2_Pm_vs_thrq_simc_norad = 0;
  TH2F *H2_Pm_vs_thrq_simc_ps = 0; //phase space
  TH2F *H2_Pm_vs_thrq_simc_ratio = 0; // for nonrad/rad ratio
  
  // define generic rootfile name to path
  TString root_file_path_rad = Form("DEUT_OUTPUT/ROOT/d2_%s_jmlfsi_rad_analyzed.root", setting.Data());
  TString root_file_path_norad = Form("DEUT_OUTPUT/ROOT/d2_%s_jmlfsi_norad_analyzed.root", setting.Data());

  
  // if file does not exist, EXIT
  if(gSystem->AccessPathName(root_file_path_rad.Data()) || gSystem->AccessPathName(root_file_path_norad.Data())) { gSystem->Exit(kTRUE);}


  // for each root file, get the desired histogram to be combined
  TFile *data_file_rad =  new TFile(root_file_path_rad.Data(), "READ");
  TFile *data_file_norad =  new TFile(root_file_path_norad.Data(), "READ");

 
  cout << Form("Read data file . . . %s",root_file_path_rad.Data()) << endl;
  
    // Retrieve histograms
  data_file_rad->cd();
  data_file_rad->GetObject("kin_plots/H_Pm_vs_thrq", H2_Pm_vs_thrq_simc_rad);  //radiative

  cout << Form("Read data file . . . %s",root_file_path_norad.Data()) << endl;
  data_file_norad->cd();
  data_file_norad->GetObject("kin_plots/H_Pm_vs_thrq", H2_Pm_vs_thrq_simc_norad); //non-radiative
  data_file_norad->GetObject("kin_plots/H_Pm_vs_thrq_ps", H2_Pm_vs_thrq_simc_ps); //non-radiative phase-space

  
  int xnbins = H2_Pm_vs_thrq_simc_norad->GetXaxis()->GetNbins();
  float xmin = H2_Pm_vs_thrq_simc_norad->GetXaxis()->GetXmin();
  float xmax = H2_Pm_vs_thrq_simc_norad->GetXaxis()->GetXmax();

  int ynbins = H2_Pm_vs_thrq_simc_norad->GetYaxis()->GetNbins();
  float ymin = H2_Pm_vs_thrq_simc_norad->GetYaxis()->GetXmin();
  float ymax = H2_Pm_vs_thrq_simc_norad->GetYaxis()->GetXmax();
  

  H2_Pm_vs_thrq_simc_ratio = new TH2F("H2_Pm_vs_thrq_simc_ratio", "SIMC Y_{norad}/Y_{rad}; Y_{norad}/Y_{rad}; #theta_{rq} [deg] ", xnbins, xmin, xmax, ynbins, ymin, ymax);
  H2_Pm_vs_thrq_simc_ratio->Divide(H2_Pm_vs_thrq_simc_norad, H2_Pm_vs_thrq_simc_rad);

      
  if(hist_type=="norad"){
    return H2_Pm_vs_thrq_simc_norad;
  }

  else if(hist_type=="rad"){
    return H2_Pm_vs_thrq_simc_rad;
  }
    
  else if(hist_type=="rad_corr_ratio"){
    return H2_Pm_vs_thrq_simc_ratio;
  }

  else if(hist_type=="phase_space"){
    return  H2_Pm_vs_thrq_simc_ps;
  }

  else{
    return dummy;
  }
  
  

}

void project2d_deut( TH2F *hist2d=0, TString setting="", Bool_t display_plots=0, Bool_t apply_radiative_corr=0 ){

  cout << "calling project2d_deut" << endl;
  //avoid display
  gROOT->SetBatch(kTRUE);
  
  /*
    Brief: Projects 2d histograms onto 1d slices in X or Y
    and plots it on a canvas subplot, and also plots relative
    errors on a separate canvas subplot

    For now is specific for deuteron, but can easily be modified for other use
   */

  // set base filename
  TString basename=Form("deut_stats_monitoring_setting_%s_", setting.Data());
  
  // Set basefilename to save projections to root file
  TString ofile="DEUT_OUTPUT/ROOT/" + basename + "output.root";
  
  TFile *fout = new TFile(ofile.Data(),"RECREATE");


  //set output names of projections to be saved to divided canvas
  TString fout_2dHist      = "DEUT_OUTPUT/PDF/" + basename + "2Dhisto.pdf";
  TString fout_projHist    = "DEUT_OUTPUT/PDF/" + basename + "projY.pdf";
  TString fout_projHistErr = "DEUT_OUTPUT/PDF/" + basename + "projY_relError.pdf";
  TString fout_projsimcRadCorr = "DEUT_OUTPUT/PDF/" + basename + "projY_simcRadCorr.pdf";
  
  // set global title/label sizes
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetLabelSize(.1, "XY");
  gStyle->SetTitleY(1.01); // offset title vertically

  // define external 2d histogrmas
  TH2F *hist2d_Pm_vs_thrq_simc_ratio = 0; //non-radiative / radiative ratio
  TH2F *hist2d_Pm_vs_thrq_simc_ps = 0;  //phase space
  TH2F *hist2d_Pm_vs_thrq_data_radcorr = 0;  //radiative corr. data
  TH2F *hist2d_Pm_vs_thrq_data_Xsec = 0;  // data cross sections
 
  cout << "Retrieve histograms get_simc_2d_histos() func. . . " << endl;

  hist2d_Pm_vs_thrq_simc_ratio = get_simc_2d_histos(setting.Data(), "rad_corr_ratio");
  hist2d_Pm_vs_thrq_simc_ps    = get_simc_2d_histos(setting.Data(), "phase_space");


  //Apply radiative corrections
  hist2d_Pm_vs_thrq_data_radcorr = (TH2F*)hist2d->Clone();
  hist2d_Pm_vs_thrq_data_radcorr->Multiply(hist2d_Pm_vs_thrq_simc_ratio);

  //scale by total charge and track inefficiencies
  
  cout << "Retrieved 2d histos . . . " << endl;
  
  // get total number of x,y bins of hist2d
  int nxbins = hist2d->GetXaxis()->GetNbins();
  int nybins = hist2d->GetYaxis()->GetNbins();


  hist2d->GetYaxis()->SetTitle("P_{m}, Missing Momentum [GeV/c]");
  hist2d->GetXaxis()->SetTitle("Recoil Angle, #theta_{rq} [deg]");

  cout << "Set title for hist2d . . . " << endl;

  hist2d_Pm_vs_thrq_simc_ratio->GetYaxis()->SetTitle("P_{m}, Missing Momentum [GeV/c]");
  hist2d_Pm_vs_thrq_simc_ratio->GetXaxis()->SetTitle("Recoil Angle, #theta_{rq} [deg]");
  hist2d_Pm_vs_thrq_simc_ratio->SetTitle("SIMC Y_{norad}/Y_{rad}");


  // plot 2d histogra
  TCanvas *c0 = new TCanvas("c0", "", 1500,1500);
  c0->Divide(2,2);

  c0->cd(1);
  hist2d->Draw("contz");
 
  
  c0->cd(2);
  gPad->Modified(); gPad->Update();
  hist2d_Pm_vs_thrq_simc_ratio->Draw("colz");

  c0->cd(3);
  gPad->Modified(); gPad->Update();
  hist2d_Pm_vs_thrq_simc_ps->Draw("colz");

  c0->cd(4);
  gPad->Modified(); gPad->Update();
  hist2d_Pm_vs_thrq_data_radcorr->Draw("colz");
  
  hist2d->Write();
  hist2d_Pm_vs_thrq_simc_ratio->Write();


  // --- define variables for calculative/plotting of relative stats. error on Pmiss (need to reset vector per projection bin) ---
  vector<double> y_val;       // this will be set to 0 (as reference)
  vector<double> y_err;       // relative error on pmiss bin counts
  vector<double> x_val;       // this is the central value of x (pmiss)
  vector<double> x_err;       // this is actually width on x-axis 
  
  double pm_counts;
  double relative_err;
  double pm_center;
  double thrq_center;
  double thrq_width;
  
  // ----------------------------

  
  // --- define canvas subplots (x,y) based on number of bins in hist2d ---
  int yc, xc;
  // if projection is along y, should make a canvas square out of nxbins
  float remainder = sqrt(nxbins) - int(sqrt(nxbins));

  if(remainder>0 && remainder<0.5){
    yc = round(sqrt(nxbins));
    xc = round(sqrt(nxbins))+1;
  }
  else{
    yc = round(sqrt(nxbins));
    xc = round(sqrt(nxbins));
  }
  
  // define canvas for projecting 2d histograms
  TCanvas *c1 = new TCanvas("c1", "Data Missing Momenta ProjY", 1400,1100);
  TCanvas *c2 = new TCanvas("c2", "Data Missing Momenta Relative Error", 1400,1100);
  TCanvas *c3 = new TCanvas("c3", "SIMC Radiative Corrections", 1400,1100);

  c1->Divide(yc, xc);
  c2->Divide(yc, xc);
  c3->Divide(yc, xc);

  
  //----------------------------


  // define 1d projection histos 
  TH1D *H_dataPm_projY = 0;
  TH1D *H_simcPm_projY_ratio = 0; // ratio of norad/rad projected
  
  cout << "About to loop over xbins of hist2d . . . " << endl;

  //loop over xbins of hist2d 
  for(int i=1; i<=nxbins; i++){

    // get xbin center value and width
    float thrq_center  = hist2d->GetXaxis()->GetBinCenter(i);
    float thrq_width   = hist2d->GetXaxis()->GetBinWidth(i); 

    //cout << "bin: " << i << ", x-val: " << thrq_center << endl;

    // project hist2d along y-axis (different bins in x)
    H_dataPm_projY       = hist2d->ProjectionY(Form("proj_Pm_thrq%.1f", thrq_center), i, i);
    H_simcPm_projY_ratio = hist2d_Pm_vs_thrq_simc_ratio->ProjectionY(Form("proj_ratio_simcPm_thrq%.1f", thrq_center), i, i);
    
    // define integrated counts on projected bin
    float counts = H_dataPm_projY->Integral();

    // reduce divisions of histos (for de-cluttering)
    H_dataPm_projY->GetYaxis()->SetNdivisions(5);
    H_dataPm_projY->GetXaxis()->SetNdivisions(10);
    H_dataPm_projY->GetXaxis()->SetLabelSize(0.1);

    H_simcPm_projY_ratio->GetYaxis()->SetNdivisions(5);
    H_simcPm_projY_ratio->GetXaxis()->SetNdivisions(10);
    H_simcPm_projY_ratio->GetXaxis()->SetLabelSize(0.1);

    //cout << "thrq_center, counts (v1) = " << thrq_center << ", " << counts << endl;
    H_dataPm_projY->SetTitle(Form("#theta_{rq} = %.1f#pm%.1f (N=%.1f)", thrq_center, thrq_width/2., counts));
    H_dataPm_projY->SetTitleSize(10);

    H_simcPm_projY_ratio->SetTitle(Form("#theta_{rq} = %.1f#pm%.1f (N=%.1f)", thrq_center, thrq_width/2., counts));
    H_simcPm_projY_ratio->SetTitleSize(10);

    //cout << Form("bin#: %d,  x-center: %.1f, counts: %.3f", i, thrq_center, counts) << endl;

    
    //-----------------------------------
    // Compute Relative Error on Counts
    //-----------------------------------
    
    // For each H_dataPm_projY hsitogram, get bin content to calculate its error and plot relative error
    y_val.clear();
    y_err.clear();
    x_val.clear();
    x_err.clear();
    
    pm_counts = 0;
    relative_err = 0;
    
    // set statistical lower (inner) limit for guidance during online data-taking
    int inner_stats = 15;  // +/- 15 %
  
    // loop over all bins of Pmiss (H_dataPm_projY histogram) 
    for(int i =1; i<=H_dataPm_projY->GetNbinsX(); i++){

      // get bin content
      pm_counts = H_dataPm_projY->GetBinContent(i);

      // get bin center
      pm_center = H_dataPm_projY->GetBinCenter(i);
      
      // calculate relative error of bin
      relative_err = (sqrt(pm_counts) / pm_counts ) * 100.; // in %

      // push values to vector
      y_val.push_back(0);
      y_err.push_back( relative_err );
      x_val.push_back(pm_center);
      x_err.push_back(0); // no need to set this width
      
      
    }

    
    // at the end, should have vector of length N for plotting
    int n=H_dataPm_projY->GetNbinsX();
    
    TGraphErrors *gr = new TGraphErrors(n, &x_val[0], &y_val[0], &x_err[0], &y_err[0]);
    gr->SetTitle(Form("proj_Pm_thrq%.1f_relErr", thrq_center));
    gr->SetMarkerColor(kBlack);
    gr->SetMarkerSize(0.);
    gr->SetMarkerStyle(21);

    

    TLine *lo_limit = new TLine( *(min_element(x_val.begin(), x_val.end())), -inner_stats, *(max_element(x_val.begin(), x_val.end())) , -inner_stats);
    TLine *up_limit = new TLine( *(min_element(x_val.begin(), x_val.end())), inner_stats, *(max_element(x_val.begin(), x_val.end())) , inner_stats);

    lo_limit->SetLineColor(kRed);
    lo_limit->SetLineStyle(1);
    lo_limit->SetLineWidth(1);

    up_limit->SetLineColor(kRed);
    up_limit->SetLineStyle(1);
    up_limit->SetLineWidth(1);
    
  
    
    //cout << "thrq_center, counts (v2) = " << thrq_center << ", " << counts << endl;
    gr->SetTitle(Form("#theta_{rq} = %.1f#pm%.1f (N=%.1f)", thrq_center, thrq_width/2., counts));
    
    // reduce divisions of histos (for de-cluttering)
    gr->GetYaxis()->SetNdivisions(5);
    gr->GetXaxis()->SetNdivisions(10);


    //---------------------------------------------------

    gPad->Modified(); gPad->Update();

      c3->cd(i);
      H_simcPm_projY_ratio->DrawClone("E0");
      H_simcPm_projY_ratio->Write();

    c1->cd(i);
    //gPad->Modified(); gPad->Update();
    H_dataPm_projY->DrawClone("histE0");
    H_dataPm_projY->Write();
        
  
    
    c2->cd(i);                                                                                                                                           
    //gPad->Modified(); gPad->Update();
    // draw to graph
    gr->Draw("AP");
    lo_limit->Draw();
    up_limit->Draw();
     
    // add legend
    if(i==1){
      auto legend = new TLegend(0.5,0.6,0.9,0.8);
      legend->AddEntry("gr",Form("#pm %d %%",inner_stats),"%d");
      legend->SetBorderSize(0);
      legend->SetTextSize(0.15);
      legend->SetTextColor(kRed);
      legend->Draw();
    }
    
    gr->Write();

    /*
  c3->cd(i);
  gPad->Modified(); gPad->Update();
    hist2d_Pm_vs_thrq_simc_ratio->DrawClone("hist");
    hist2d_Pm_vs_thrq_simc_ratio->Write();
    */
    
    
  } // end loop over 2D xbins [th_rq]

  
  // save canvas
  gStyle->SetOptStat(0);
  c0->SaveAs( fout_2dHist.Data()      );
  c1->SaveAs( fout_projHist.Data()    );
  c2->SaveAs( fout_projHistErr.Data() );
  c3->SaveAs( fout_projsimcRadCorr.Data() );

  if(display_plots){
    
    //open plots with evince or any other viewer
    /*
    gSystem->Exec(Form("evince %s", fout_2dHist.Data() ));
    gSystem->Exec(Form("evince %s", fout_projHist.Data() ));
    gSystem->Exec(Form("evince %s", fout_projHistErr.Data() )); 
    gSystem->Exec(Form("evince %s",  fout_projsimcRadCorr.Data() ));
    */

    gSystem->Exec(Form("open %s", fout_2dHist.Data() ));
    gSystem->Exec(Form("open %s", fout_projHist.Data() ));
    gSystem->Exec(Form("open %s", fout_projHistErr.Data() )); 
    gSystem->Exec(Form("open %s",  fout_projsimcRadCorr.Data() ));
    
  }
  
}


void project2d_online_develop() {

  int run_min=20871;
  int run_max=20871;
  TString pm_setting="pm120";
  /*
  cout << "" << endl;
  cout << "----------------------------------------------------------------------------" << endl;
  cout << "" << endl;
  cout << "function: project2d_online(int run_min, int run_max, TString pm_setting)" << endl;
  cout << "" << endl;
  cout << "Brief: This function plots the following: \n"
  "(1) the combined 2d histogram Pmiss vs. th_rq \n"
  "(2) projection of Pmiss (y-axis) in th_rq bins (x-axis)   \n"
  "(3) relative statistical errors (sqrt[N]/[N]) of projected bins\n (for monitoring statistical goals per thrq bin)" << endl;
  cout << "" << endl;
  cout <<  "**NOTE** : if a run range is over a single or multiple settings, \n "
    "name it accordintly e.g. pm_120 (if single), pm_120_580 (multiple settings) etc. \n "
    "This will be used to name the output .pdf \n" << endl;
  cout << "" << endl;
  cout << "----------------------------------------------------------------------------" << endl;
  cout << "" << endl;
  cout << "Please enter minimum run in range (e.g., 3288): ";
  cin >> run_min;
  cout << "\n Please enter maximum run in range (e.g., 3377): ";
  cin >> run_max;
  cout << "\n Please enter Pmiss setting (e.g. pm_120, pm_120_580): ";
  cin >> pm_setting;
  cout << "" << endl;
  cout << "----------------------------------------------------------------------------" << endl;
  */
  // Brief: this function calls two functions described below. 

  
  // funct1: retrieve the combined 2d histos on a given run range. The root file is predefined in the function,
  // so  all the user needs is to input the path where the histogram is on the .root file
  TH2F * myhist2d = combine_2dhistos(run_min, run_max, "randSub_plots/H_Pm_vs_thrq_rand_sub");

  // func2: projectes the 2d histo (slices of x-bins along y-axis) onto 1d bins, the pm_setting is just for histogram naming purposes
  // and should be consistent with the histogram range chosen
  project2d_deut( myhist2d, pm_setting, true, true ); // the bool flag is to display the plots (otherwise, they will be saved)

}

//#endif
