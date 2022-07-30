// void pretty_TH2D(TH2D *h2);
// void pretty_TH1D(TH1D *h2);
void plotOptics_modified(int runNUM = -1, int evtNUM = 0, TString ana_type = "") {

  // TFile *f = new
  // TFile(Form("/volatile/hallc/spring17/holly/shmsoptics/shms_replay_production_all_%i_-1.root",runNUM));

  //Dien's testing
  //TString ifname_pattern =
  //    Form("./cafe_replay_%s_%d_%d.root", ana_type.Data(), runNUM, evtNUM);

  // C.Y. set input ROOTfile Name Pattern for CaFe
  TString ifname_pattern=Form("ROOTfiles/%s/cafe_replay_%s_%d_%d.root",
			      ana_type.Data(), ana_type.Data(), runNUM, evtNUM);

  TFile *f = new TFile(ifname_pattern.Data());

  TTree *tt = (TTree *)f->Get("T");

  // here's the cut
  TCut cut =
      "P.gtr.dp<20&&P.cal.etracknorm>0.8&&P.ngcer.npeSum>5";

  TCut z_pos = "P.react.z > 5 && P.react.z <10";
  TCut z_neg = "P.react.z > -11 && P.react.z < -6 ";
  TCut delta_0 = "P.gtr.dp > 0";  // need to set this to > 0% for cafe

  TCut cutCentral =
      "abs(P.gtr.x+P.gtr.th*253.0)<1&&abs((-0.019*P.gtr.dp+0.00019*P.gtr.dp*P."
      "gtr.dp+(138.0+75.0)*P.gtr.ph+P.gtr.y) + "
      "40.0*(-0.00052*P.gtr.dp+0.0000052*pow(P.gtr.dp,2)+P.gtr.ph))<1";

  // make the output file
  TCanvas *canvas = new TCanvas("canvas", "plots", 800, 600);
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(0);
  canvas->SetFrameFillColor(0);
  canvas->SetFrameBorderMode(0);

  gROOT->SetBatch(false);
  gStyle->SetOptStat(0);
  canvas->SetGridx();
  canvas->SetGridy();

  std::string pdf_file_name =
      Form("CAFE_OUTPUT/PDF/cafe_optics_%s_%i_%i_sieve0.pdf", ana_type.Data(),
           runNUM, evtNUM);

  TFile *fout =
      new TFile(Form("CAFE_OUTPUT/ROOT/cafe_optics_%s_%i_%i_sieve0.root",
                     ana_type.Data(), runNUM, evtNUM),
                "RECREATE");

  // Define histogram for plots without the z-vertex cut
  TH1F *h_z = new TH1F("h_z", ";P.react.z [cm]", 200, -15, 15);
  TH2F *h2_ypVy =
      new TH2F("h2_ypVy", ";yTar [cm];ypTar", 100, -5, 5, 100, -0.05, 0.05);
  TH2F *h2_ypVz = new TH2F("h2_ypVz", ";zVertex [cm];ypTar", 100, -14, 14, 100,
                           -0.05, 0.05);
  TH2F *h2_dpVz =
      new TH2F("h2_dpVz", ";zVertex [cm];delta", 100, -14, 14, 100, -15, 20);
  TH2F *h2_ypVxp = new TH2F("h2_ypVxp", "; xpTar; ypTar", 100, -0.05, 0.05, 120,
                            -0.06, 0.06);
  TH2F *h2_yfpVxfp =
      new TH2F("h2_yfpVxfp", ";xfp [cm];yfp [cm]", 100, 0, 8, 100, -10, 10);

  TH2F *h2_sieve = new TH2F("h2_sieve", ";ySieve [cm];xSieve[cm]", 180, -9.0,
                            9.0, 260, -13.0, 13.0);
  TH2F *h2_xpVd =
      new TH2F("h2_xpVd", ";delta;xpfp", 100, -10, 20, 100, -0.15, 0.15);

  // Define histogram for plots for pos foil
  TH2F *h2_sieve_pos =
      new TH2F("h2_sieve_pos",
               "Sieve for Positive Vz, delta >0;ysieve [cm]; xSieve [cm]", 180,
               -9, 9, 260, -13, 13);

  // Define histogram for plots for neg foil
  TH2F *h2_sieve_neg =
      new TH2F("h2_sieve_neg",
               "Sieve for Negative Vz, delta >0 ;ysieve [cm]; xSieve [cm]", 180,
               -9, 9, 260, -13, 13);

  // Define histogram for plots with central hole only
  TH1F *h_z_c =
      new TH1F("h_z_c", "central sieve hole;P.react.z [cm]", 100, -14, 14);
  TH2F *h2_ypVy_c = new TH2F("h2_ypVy_c", "central sieve hole;yTar [cm];ypTar",
                             100, -5, 5, 100, -0.05, 0.05);
  TH2F *h2_yfpVxfp_c =
      new TH2F("h2_yfpVxfp_c", "central sieve hole;yfp [cm];yfp [cm]", 100, 0,
               8, 100, -10, 10);
  TH2F *h2_dpVz_c =
      new TH2F("h2_dpVz_c", "central sieve hole;zVertex [cm];delta", 100, -14,
               14, 100, -10, 20);
  TH2F *h2_ypVz_c =
      new TH2F("h2_ypVz_c", "central sieve hole;zVertex [cm];ypTar", 100, -14,
               14, 100, -0.05, 0.05);
  TH2F *h2_sieve_c =
      new TH2F("h2_sieve_c", "central sieve hole;ySieve [cm];xSieve[cm]", 200,
               -7.0, 7.0, 200, -12.0, 12.0);
  TH2F *h2_xpVd_c = new TH2F("h2_xpVd_c", "central sieve hole;delta;xpfp", 100,
                             -15, 20, 100, -0.15, 0.15);

  TH2F *h2_ypVzSlice[8];

  // plot this stuff
  tt->Draw("P.gtr.ph:P.react.z>>h2_ypVz", cut);
  tt->Draw("P.dc.y_fp:P.dc.x_fp>>h2_yfpVxfp", cut);
  tt->Draw("P.gtr.ph:P.gtr.y>>h2_ypVy", cut);
  tt->Draw("P.react.z>>h_z", cut);
  tt->Draw("P.gtr.dp:P.react.z>>h2_dpVz", cut);
  tt->Draw("P.gtr.ph:P.gtr.th>>h2_ypVxp", cut);
  tt->Draw(
      "P.gtr.x+P.gtr.th*253.0:(-0.019*P.gtr.dp+0.00019*P.gtr.dp*P.gtr.dp+(138."
      "0+75.0)*P.gtr.ph+P.gtr.y) + "
      "40.0*(-0.00052*P.gtr.dp+0.0000052*pow(P.gtr.dp,2)+P.gtr.ph)>>h2_sieve",
      cut);
  tt->Draw("P.dc.xp_fp:P.gtr.dp>>h2_xpVd", cut);

  // Plotting for central hole
  tt->Draw("P.gtr.ph:P.react.z>>h2_ypVz_c", cut && cutCentral);
  tt->Draw("P.dc.y_fp:P.dc.x_fp>>h2_yfpVxfp_c", cut && cutCentral);
  tt->Draw("P.gtr.ph:P.gtr.y>>h2_ypVy_c", cut && cutCentral);
  tt->Draw("P.react.z>>h_z_c", cut && cutCentral);
  tt->Draw("P.gtr.dp:P.react.z>>h2_dpVz_c", cut && cutCentral);
  tt->Draw(
      "P.gtr.x+P.gtr.th*253.0:(-0.019*P.gtr.dp+0.00019*P.gtr.dp*P.gtr.dp+(138."
      "0+75.0)*P.gtr.ph+P.gtr.y) + "
      "40.0*(-0.00052*P.gtr.dp+0.0000052*pow(P.gtr.dp,2)+P.gtr.ph)>>h2_sieve_c",
      cut && cutCentral);
  tt->Draw("P.dc.xp_fp:P.gtr.dp>>h2_xpVd_c", cut && cutCentral);

  for (int ii = 0; ii < 8; ii++) {
    h2_ypVzSlice[ii] =
        new TH2F(Form("h2_ypVzSlice_%i", ii),
                 Form("xfp=%i cm +/- 0.5cm;zVertex [cm]; ypTar", ii + 1), 100,
                 -15, 15, 100, -0.05, 0.05);
    TCut slice = Form("abs(P.dc.x_fp - (%i+1))<0.5", ii);
    tt->Draw(Form("P.gtr.ph:P.react.z>>h2_ypVzSlice_%i", ii), cut && slice);
  }

  // Plot the 2D sieve for pos foil
  tt->Draw(
      "P.gtr.x+P.gtr.th*253.0:(-0.019*P.gtr.dp+0.00019*P.gtr.dp*P.gtr.dp+(138."
      "0+75.0)*P.gtr.ph+P.gtr.y) + "
      "40.0*(-0.00052*P.gtr.dp+0.0000052*pow(P.gtr.dp,2)+P.gtr.ph)>>h2_sieve_"
      "pos",
      cut && z_pos && delta_0, "colz");

  // plot the 2D sieve for neg foil
  tt->Draw(
      "P.gtr.x+P.gtr.th*253.0:(-0.019*P.gtr.dp+0.00019*P.gtr.dp*P.gtr.dp+(138."
      "0+75.0)*P.gtr.ph+P.gtr.y) + "
      "40.0*(-0.00052*P.gtr.dp+0.0000052*pow(P.gtr.dp,2)+P.gtr.ph)>>h2_sieve_"
      "neg",
      cut && z_neg && delta_0, "colz");

  // save plots
  canvas->Update();

  h_z->SetLineWidth(2);
  h_z->Draw();
  canvas->Print((pdf_file_name + "(").c_str());

  h2_ypVz->Draw("colz");
  canvas->Print((pdf_file_name + "(").c_str());

  h2_ypVy->Draw("colz");
  canvas->Print((pdf_file_name + "(").c_str());

  h2_dpVz->Draw("colz");
  canvas->Print((pdf_file_name + "(").c_str());

  h2_ypVxp->Draw("colz");
  canvas->Print((pdf_file_name + "(").c_str());

  h2_yfpVxfp->Draw("colz");
  canvas->Print((pdf_file_name + "(").c_str());

  h2_sieve->Draw("colz");
  canvas->Print((pdf_file_name + "(").c_str());

  h2_xpVd->Draw("colz");
  canvas->Print((pdf_file_name + "(").c_str());

  // adding the sieve plots for pos and neg file here
  h2_sieve_pos->Draw("colz");
  canvas->Print((pdf_file_name + "(").c_str());

  h2_sieve_neg->Draw("colz");
  canvas->Print((pdf_file_name + ")").c_str());

  // Adding the central hole
  // h_z_c->SetLineWidth(2);
  // h_z_c->Draw();
  // canvas->Print((pdf_file_name + "(").c_str());

  // h2_ypVz_c->Draw("colz");
  // canvas->Print((pdf_file_name + "(").c_str());

  // h2_dpVz_c->Draw("colz");
  // canvas->Print((pdf_file_name + "(").c_str());

  // h2_ypVy_c->Draw("colz");
  // canvas->Print((pdf_file_name + "(").c_str());

  // h2_yfpVxfp_c->Draw("colz");
  // canvas->Print((pdf_file_name + "(").c_str());

  // h2_sieve_c->Draw("colz");
  // canvas->Print((pdf_file_name + "(").c_str());

  // h2_xpVd_c->Draw("colz");
  // canvas->Print((pdf_file_name + "(").c_str());

  // for (int jj = 0; jj < 8; jj++) {
  //   h2_ypVzSlice[jj]->Draw("colz");
  //  canvas->Print((pdf_file_name + "(").c_str());
  //}
  // h2_xpVd_c->Draw("colz");
  // canvas->Print((pdf_file_name + ")").c_str());

  fout->Write();
  fout->Close();

  // gSystem->Exec(Form("evince %s", pdf_file_name.c_str()));
}
