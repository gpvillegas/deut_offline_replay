#include <vector>


void project2d( TH2F *hist2d=0 ){
//void project2d( ){

  /*
    Brief: Projects 2d histograms onto 1d slices in X or Y
    and plots it on a canvas subplot, and also plots relative
    errors on a separate canvas subplot
   */
  
  // set global title/label sizes
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetLabelSize(.1, "XY");
  gStyle->SetTitleY(1.01); // offset title vertically

  //TH2F *hist2d = new TH2F("hist2d", "", 30,-15,15, 30, -15,15);


  // get total number of x,y bins of hist2d
  int nxbins = hist2d->GetXaxis()->GetNbins();
  int nybins = hist2d->GetYaxis()->GetNbins();

  
  
  //TF2 *xyg = new TF2("xyg","xygaus",-10,10,-10,10); xyg->SetParameters(1,5,2,5,2);
  //hist2d->FillRandom("xyg");

  // plot 2d histogram
  TCanvas *c0 = new TCanvas("c0", "", 1200,800);
  c0->cd();
  /*
  hist2d->SetTitle("P_{miss} vs. #theta_{rq}");
  hist2d->SetTitleSize(0.03);
  hist2d->GetYaxis()->SetLabelSize(0.05);
  hist2d->GetXaxis()->SetLabelSize(0.05);
  hist2d->GetYaxis()->SetTitle("P_{m} [GeV/c]");
  hist2d->GetXaxis()->SetTitle("#theta_{rq} [deg]");
  */
  
  hist2d->Draw("colz");

  
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
  yc = round(sqrt(nxbins));
  xc = yc+1;

  TCanvas *c1 = new TCanvas("c1", "", 1400,1100);
  TCanvas *c2 = new TCanvas("c2", "", 1400,1100);
  
  c1->Divide(yc, xc);
  c2->Divide(yc, xc);
  //----------------------------

  
  //loop over xbins of hist2d 
  for(int i=1; i<=nxbins; i++){

    // get xbin center value and width
    float thrq_center  = hist2d->GetXaxis()->GetBinCenter(i);
    float thrq_width   = hist2d->GetXaxis()->GetBinWidth(i); 

    //cout << "bin: " << i << ", x-val: " << thrq_center << endl;

    // project hist2d along y-axis (different bins in x)
    TH1D *h1 = hist2d->ProjectionY("projy", i, i);

    // define integrated counts on projected bin
    float counts = h1->Integral();

    // change to canvas subplot of ith bin projection
    c1->cd(i);

    // reduce divisions of histos (for de-cluttering)
    h1->GetYaxis()->SetNdivisions(5);
    h1->GetXaxis()->SetNdivisions(10);
    h1->Draw("histE0");
    //cout << "thrq_center, counts (v1) = " << thrq_center << ", " << counts << endl;
    h1->SetTitle(Form("#theta_{rq} = %.1f#pm%.1f (N=%.1f)", thrq_center, thrq_width/2., counts));
    h1->SetTitleSize(10);
   
    c1->Update();
    
    //cout << Form("bin#: %d,  x-center: %.1f, counts: %.3f", i, xbin_center, counts) << endl;


    //-----------------------------------
    // Compute Relative Error on Counts
    //-----------------------------------
    
    // For each h1 hsitogram, get bin content to calculate its error and plot relative error
    y_val.clear();
    y_err.clear();
    x_val.clear();
    x_err.clear();
    
    pm_counts = 0;
    relative_err = 0;
    
    // loop over all bins of Pmiss (h1 histogram) 
    for(int i =1; i<=h1->GetNbinsX(); i++){

      // get bin content
      pm_counts = h1->GetBinContent(i);

      // get bin center
      pm_center = h1->GetBinCenter(i);
      
      // calculate relative error of bin
      relative_err = (sqrt(pm_counts) / pm_counts ) * 100.; // in %

      // push values to vector
      y_val.push_back(0);
      y_err.push_back( relative_err );
      x_val.push_back(pm_center);
      x_err.push_back(0); // no need to set this width
      
      
    }

    // at the end, should have vector of length N for plotting
    int n=h1->GetNbinsX();
    
    TGraphErrors *gr = new TGraphErrors(n, &x_val[0], &y_val[0], &x_err[0], &y_err[0]);
 
    gr->SetMarkerColor(kBlack);
    gr->SetMarkerSize(0.);
    gr->SetMarkerStyle(21);

    c2->cd(i);

    TLine *lo_limit = new TLine( *(min_element(x_val.begin(), x_val.end())), -15., *(max_element(x_val.begin(), x_val.end())) , -15.);
    TLine *up_limit = new TLine( *(min_element(x_val.begin(), x_val.end())), 15., *(max_element(x_val.begin(), x_val.end())) , 15.);

    lo_limit->SetLineColor(kBlack);
    lo_limit->SetLineStyle(1);
    lo_limit->SetLineWidth(1);

    up_limit->SetLineColor(kBlack);
    up_limit->SetLineStyle(1);
    up_limit->SetLineWidth(1);
    
  
    
    cout << "thrq_center, counts (v2) = " << thrq_center << ", " << counts << endl;
    gr->SetTitle(Form("#theta_{rq} = %.1f#pm%.1f (N=%.1f)", thrq_center, thrq_width/2., counts));
    
    // reduce divisions of histos (for de-cluttering)
    gr->GetYaxis()->SetNdivisions(5);
    gr->GetXaxis()->SetNdivisions(10);
    
    // draw to graph
    gr->Draw("AP");
    lo_limit->Draw();
    up_limit->Draw("same");
    c2->Update();

    
      
  } // end loop over 2D xbins [th_rq]
  
  
    
    
}
