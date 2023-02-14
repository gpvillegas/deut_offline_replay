#include "TCanvas.h"
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TLine.h"
#include "TProfile.h"
#include "TH2F.h"
#include <iostream>

#ifndef CHECK_CALIB
#define CHECK_CALIB


//=====================================================
//==============DEFINE SOME CONSTANTS==================
//=====================================================

//Define some detectors planes and sides
static const Int_t hod_PLANES = 4;
static const Int_t cal_PLANES = 4;
static const Int_t dc_PLANES = 12;
static const Int_t SIDES = 2;

static const string hod_pl_names[hod_PLANES] = {"1x", "1y", "2x", "2y"};
static const string cal_pl_names[cal_PLANES] = {"1pr", "2ta", "3ta", "4ta"};

static const string hdc_pl_names[dc_PLANES] = {"1u1", "1u2", "1x1", "1x2", "1v2", "1v1", "2v1", "2v2", "2x2", "2x1", "2u2", "2u1"};
static const string pdc_pl_names[dc_PLANES] = {"1u1", "1u2", "1x1", "1x2", "1v1", "1v2", "2v2", "2v1", "2x2", "2x1", "2u2", "2u1"};
 
static const Int_t hnwires[dc_PLANES] = {96, 96, 102, 102, 96, 96, 96, 96, 102, 102, 96, 96};
static const Int_t pnwires[dc_PLANES] = {107, 107, 79, 79, 107, 107, 107, 107, 79, 79, 107, 107};

static const string side_names[SIDES] = {"GoodPos", "GoodNeg"};
static const string cal_side_names[SIDES] = {"goodPos", "goodNeg"};

static const string nsign[SIDES] = {"+", "-"};

static const Int_t hmaxPMT[hod_PLANES] = {16, 10, 16, 10};
static const Int_t pmaxPMT[hod_PLANES] = {13, 13, 14, 21};


//=====================================================



//=============================================
//===========DEFINE HISTOGRAMS=================
//=============================================

//Set Histogram Binning

// HMS Focal Plane
Double_t hxfp_nbins,   hxfp_xmin, hxfp_xmax;
Double_t hyfp_nbins,   hyfp_xmin, hyfp_xmax;   
Double_t hdelta_nbins, hdelta_xmin, hdelta_xmax;   

// sHMS Focal Plane 
Double_t pxfp_nbins,   pxfp_xmin,   pxfp_xmax; 
Double_t pyfp_nbins,   pyfp_xmin,   pyfp_xmax;  
Double_t pdelta_nbins, pdelta_xmin, pdelta_xmax;   

//HMS DC
Double_t hdcTime_nbins, hdcTime_xmin, hdcTime_xmax;
Double_t hdcDist_nbins, hdcDist_xmin, hdcDist_xmax;
Double_t hdcRes_nbins, hdcRes_xmin, hdcRes_xmax;

//HMS HODO
Double_t hhodBeta_nbins, hhodBeta_xmin, hhodBeta_xmax;
Double_t hhodXtrk_nbins, hhodXtrk_xmin, hhodXtrk_xmax;  
Double_t hhodYtrk_nbins, hhodYtrk_xmin, hhodYtrk_xmax;  
//HMS CALO
Double_t hcalEtrkNorm_nbins, hcalEtrkNorm_xmin, hcalEtrkNorm_xmax;
Double_t hcalEtot_nbins, hcalEtot_xmin, hcalEtot_xmax;
Double_t hcalXtrk_nbins, hcalXtrk_xmin, hcalXtrk_xmax;  
Double_t hcalYtrk_nbins, hcalYtrk_xmin, hcalYtrk_xmax; 
//HMS CER
Double_t hcer_nbins, hcer_xmin, hcer_xmax;

//SHMS DC
Double_t pdcTime_nbins, pdcTime_xmin, pdcTime_xmax;
Double_t pdcDist_nbins, pdcDist_xmin, pdcDist_xmax;
Double_t pdcRes_nbins, pdcRes_xmin, pdcRes_xmax;
//SHMS HODO
Double_t phodBeta_nbins, phodBeta_xmin, phodBeta_xmax;
Double_t phodXtrk_nbins, phodXtrk_xmin, phodXtrk_xmax;  
Double_t phodYtrk_nbins, phodYtrk_xmin, phodYtrk_xmax;  
//SHMS CALO
Double_t pcalEtrkNorm_nbins, pcalEtrkNorm_xmin, pcalEtrkNorm_xmax;
Double_t pcalEtot_nbins, pcalEtot_xmin, pcalEtot_xmax;
Double_t pcalXtrk_nbins, pcalXtrk_xmin, pcalXtrk_xmax;  
Double_t pcalYtrk_nbins, pcalYtrk_xmin, pcalYtrk_xmax; 
//SHMS HGCER/NGCER
Double_t phgcer_nbins, phgcer_xmin, phgcer_xmax;
Double_t pngcer_nbins, pngcer_xmin, pngcer_xmax;


//=========================
//====DEFINE CANVAS========
//=========================

//====HMS=======                            //====SHMS=====
                                             				       
//==HODOSCOPES===		             //==HODOSCOPES===		                     
TCanvas *hhodCanv;		       	     TCanvas *phodCanv;		       
TCanvas *hhodCanv2D;		       	     TCanvas *phodCanv2D;		       
TCanvas *hhodProfCanv;		       	     TCanvas *phodProfCanv;		       
				       	     				       
//==DRIFT CHAMBERS===		       	     //==DRIFT CHAMBERS===		       
TCanvas *hdcTimeCanv;		       	     TCanvas *pdcTimeCanv;		       
TCanvas *hdcDistCanv;		       	     TCanvas *pdcDistCanv;		       
TCanvas *hdcResCanv;		       	     TCanvas *pdcResCanv;		       
TCanvas *hdcResGraphCanv;	       	     TCanvas *pdcResGraphCanv;	       
				       	     				       
//Canvas to Draw 2D variales 	       	     //Canvas to Draw 2D variales 	       
//vs. wire num and Profile Histos      	     //vs. wire num and Profile Histos      
TCanvas *hdcResCanv2D;		       	     TCanvas *pdcResCanv2D;		       
TCanvas *hdcTimeCanv2D;		       	     TCanvas *pdcTimeCanv2D;		       
TCanvas *hdcDistCanv2D;		       	     TCanvas *pdcDistCanv2D;		       
TCanvas *hdcResCanvProf;	       	     TCanvas *pdcResCanvProf;	       
				       	     				       
//==CALORIMETER===		       	     //==CALORIMETER===		       
TCanvas *hcalCanv;		       	     TCanvas *pcalCanv;		       
				       	     				       
//==CHERENKOV===		       	     //==CHERENKOV===		       
TCanvas *hcerCanv;                     	     TCanvas *phgcerCanv; TCanvas *pngcerCanv;                     
                                             TCanvas *phgcerSumCanv; TCanvas *pngcerSumCanv;   

//===================
//Define Histograms
//==================

//=====HMS=======                                    //=====SHMS======

//===DRIFT CHAMBERS====			             //===DRIFT CHAMBERS====			       
TH1F *H_hdcDist[dc_PLANES];		      	     TH1F *H_pdcDist[dc_PLANES];		      
TH1F *H_hdcTime[dc_PLANES];		      	     TH1F *H_pdcTime[dc_PLANES];		      
TH1F *H_hdcRes[dc_PLANES];		      	     TH1F *H_pdcRes[dc_PLANES];		      
					      	     					      
TH2F *H_hres_vs_wire[dc_PLANES];       	      	     TH2F *H_pres_vs_wire[dc_PLANES];		      
TH2F *H_htime_vs_wire[dc_PLANES];	      	     TH2F *H_ptime_vs_wire[dc_PLANES];	      
TH2F *H_hdist_vs_wire[dc_PLANES];	      	     TH2F *H_pdist_vs_wire[dc_PLANES];	      
					      	     					      
TProfile *hdcResProf[dc_PLANES];	       	     TProfile *pdcResProf[dc_PLANES];		      
					      	     					      
//==CALORIMETER HISTOGRAMS===		      	     //==CALORIMETER HISTOGRAMS===		      
TH1F *H_hcalEtrkNorm;			      	     TH1F *H_pcalEtrkNorm;			      
TH1F *H_hcalEtot;			      	     TH1F *H_pcalEtot;			      
TH2F *H_hcalEtrkNorm_vs_xtrk;		      	     TH2F *H_pcalEtrkNorm_vs_xtrk;		      
TH2F *H_hcalEtrkNorm_vs_ytrk;		      	     TH2F *H_pcalEtrkNorm_vs_ytrk;		      

TH2F *H_hcalEtrkNorm_vs_xfp;                         TH2F *H_pcalEtrkNorm_vs_xfp;
TH2F *H_hcalEtrkNorm_vs_yfp;                         TH2F *H_pcalEtrkNorm_vs_yfp;
TH2F *H_hcalEtrkNorm_vs_delta;                       TH2F *H_pcalEtrkNorm_vs_delta;

					      	     					      
//===HODOSCOPES===			      	     //===HODOSCOPES===			      
TH1F *H_hbeta_peak;                                  TH1F *H_pbeta_peak;
TH1F *H_hhodBeta;			      	     TH1F *H_phodBeta;			      
TH1F *H_hhodBetaNoTrk;			      	     TH1F *H_phodBetaNoTrk;			      
TH2F *H_hhodBeta_v_Xtrk[hod_PLANES];	      	     TH2F *H_phodBeta_v_Xtrk[hod_PLANES];	      
TH2F *H_hhodBeta_v_Ytrk[hod_PLANES];	      	     TH2F *H_phodBeta_v_Ytrk[hod_PLANES];	      

TH2F *H_hhodBeta_vs_xfp;                              TH2F *H_phodBeta_vs_xfp; 
TH2F *H_hhodBeta_vs_yfp;                              TH2F *H_phodBeta_vs_yfp; 
TH2F *H_hhodBeta_vs_delta;                            TH2F *H_phodBeta_vs_delta;

					      	     					      
//Profile Histos of beta vs. xtrk (or ytrk)   	     //Profile Histos of beta vs. xtrk (or ytrk)   
TProfile *hhod_xProfX[hod_PLANES];	      	     TProfile *phod_xProfX[hod_PLANES];	      
TProfile *hhod_yProfX[hod_PLANES];	      	     TProfile *phod_yProfX[hod_PLANES];	      
					      	     					      
//===CHERENKOV===			      	     //===CHERENKOV===			      
TH1F *H_hcerNpe[2];			      	     TH1F *H_phgcerNpe[4];   TH1F *H_pngcerNpe[4];	    
TH1F *H_hcerNpeSum;                           	     TH1F *H_phgcerNpeSum;   TH1F *H_pngcerNpeSum;          


//========================
//Define TTree Leaf Names 
//========================

//---Names---
TString base;
TString nhdc_wire;			          TString npdc_wire;			       
TString nhdc_time;			       	  TString npdc_time;			       
TString nhdc_ndata;			       	  TString npdc_ndata;			       
TString nhdc_dist;			       	  TString npdc_dist;			       
TString nhdc_nhit;			       	  TString npdc_nhit;			       
TString nhdc_res;			       	  TString npdc_res;			       
					       	  					       
TString nhhod_beta;			       	  TString nphod_beta;			       
TString nhhod_beta_notrk;		       	  TString nphod_beta_notrk;		       
TString nhhod_xtrack;			       	  TString nphod_xtrack;			       
TString nhhod_ytrack;			       	  TString nphod_ytrack;			       
					       	  					       
TString nhcer_npesum;			       	  TString nphgcer_npesum;   TString npngcer_npesum;			       
TString nhcer_npe;			       	  TString nphgcer_npe;	    TString npngcer_npe;		       
					       	  					       
TString nhcal_etot;	     		       	  TString npcal_etot;	     		       
TString nhcal_etrknorm; 		       	  TString npcal_etrknorm; 		       
TString nhcal_xtrack;	     		       	  TString npcal_xtrack;	     		       
TString nhcal_ytrack;                          	  TString npcal_ytrack;                          

TString nhxfp;       TString npxfp;
TString nhyfp;       TString npyfp;
TString nhdelta;     TString npdelta; 


//--Variables
Double_t hdc_wire[dc_PLANES][1000];	     Double_t pdc_wire[dc_PLANES][1000];	   
Double_t hdc_time[dc_PLANES][1000];	     Double_t pdc_time[dc_PLANES][1000];	   
Int_t hdc_ndata[dc_PLANES];		     Int_t pdc_ndata[dc_PLANES];		   
Double_t hdc_dist[dc_PLANES][1000];	     Double_t pdc_dist[dc_PLANES][1000];	   
Double_t hdc_res[dc_PLANES]; 		     Double_t pdc_res[dc_PLANES]; 		   
Double_t hdc_nhit[dc_PLANES];		     Double_t pdc_nhit[dc_PLANES];		   
					     					   
Double_t hhod_beta;			     Double_t phod_beta;			   
Double_t hhod_beta_notrk;		     Double_t phod_beta_notrk;		   
Double_t hhod_xtrack[hod_PLANES];	     Double_t phod_xtrack[hod_PLANES];	   
Double_t hhod_ytrack[hod_PLANES];	     Double_t phod_ytrack[hod_PLANES];	   
					     					   
Double_t hcer_npe[2];			     Double_t phgcer_npe[2];    Double_t pngcer_npe[2];			   
Double_t hcer_npesum;			     Double_t phgcer_npesum;	Double_t pngcer_npesum;		   
					     					   
Double_t hcal_etot;			     Double_t pcal_etot;			   
Double_t hcal_etrknorm;			     Double_t pcal_etrknorm;			   
Double_t hcal_xtrack;			     Double_t pcal_xtrack;			   
Double_t hcal_ytrack;                        Double_t pcal_ytrack;                      

Double_t hxfp;      Double_t pxfp;
Double_t hyfp;      Double_t pyfp;  
Double_t hdelta;    Double_t pdelta;  
  
//Define Mean/Sigma to be used for residuals
Double_t mean[dc_PLANES];
Double_t mean_err[dc_PLANES];

Double_t sigma[dc_PLANES];
Double_t sigma_err[dc_PLANES];

Double_t x[dc_PLANES];
Double_t x_err[dc_PLANES];


//HMS Detector Leaf Names                                                                                          


//-------------------------------------------------------

//--------Define Histograms for reference time Cuts------

//HMS Histograms
TH1F *H_hodo_Tref;
TH1F *H_DC_Tref[4];
TH1F *H_FADC_Tref;

//HMS Multiplicity CUT Histograms (Apply mult == 1 cut for HMS ref. times, as there was only 1 ref. time, 3/4 for protons)
TH1F *H_hodo_Tref_CUT;
TH1F *H_DC_Tref_CUT[4];
TH1F *H_FADC_Tref_CUT;

//SHMS Histograms
TH1F *P_hodo_Tref1;
TH1F *P_hodo_Tref2;
TH1F *P_DC_Tref[10];
TH1F *P_FADC_Tref;

//SHMS Multiplicity CUT Histograms (Applt mult == 3 cut for SHMS ref. time, as there were 3 ref. time, 3/4, pEL_REAL, pEL_CLEAN)
TH1F *P_hodo_Tref1_CUT;
TH1F *P_hodo_Tref2_CUT;
TH1F *P_DC_Tref_CUT[10];
TH1F *P_FADC_Tref_CUT;

//-------Define Histograms for Time Window Cuts--------

//HMS Histograms
TH1F *H_hod_TdcAdcTimeDiff[hod_PLANES][SIDES][16];
TH1F *H_cal_TdcAdcTimeDiff[cal_PLANES][SIDES][13];
TH1F *H_dc_rawTDC[dc_PLANES];
TH1F *H_cer_TdcAdcTimeDiff[2];

//HMS Histograms with Multiplicity == 1 CUT
TH1F *H_hod_TdcAdcTimeDiff_CUT[hod_PLANES][SIDES][16];
TH1F *H_cal_TdcAdcTimeDiff_CUT[cal_PLANES][SIDES][13];
TH1F *H_dc_rawTDC_CUT[dc_PLANES];
TH1F *H_cer_TdcAdcTimeDiff_CUT[2];

//SHMS Histograms
TH1F *P_hod_TdcAdcTimeDiff[hod_PLANES][SIDES][21];
TH1F *P_cal_TdcAdcTimeDiff[224];  //fly's eye (224 pmt-channels)
TH1F *P_prSh_TdcAdcTimeDiff[SIDES][14];
TH1F *P_dc_rawTDC[dc_PLANES];
TH1F *P_hgcer_TdcAdcTimeDiff[4];
TH1F *P_ngcer_TdcAdcTimeDiff[4];

//SHMS Histograms with Multiplicity == 3 CUT
TH1F *P_hod_TdcAdcTimeDiff_CUT[hod_PLANES][SIDES][21];
TH1F *P_cal_TdcAdcTimeDiff_CUT[224];  
TH1F *P_prSh_TdcAdcTimeDiff_CUT[SIDES][14];
TH1F *P_dc_rawTDC_CUT[dc_PLANES];
TH1F *P_hgcer_TdcAdcTimeDiff_CUT[4];
TH1F *P_ngcer_TdcAdcTimeDiff_CUT[4];

//TRG DETECTOR Histograms
TH1F *pTrig1_ROC1_rawTdcTime;
TH1F *pTrig1_ROC2_rawTdcTime;
TH1F *pTrig4_ROC1_rawTdcTime;
TH1F *pTrig4_ROC2_rawTdcTime;

//=========================
//====DEFINE CANVAS========
//=========================

//Define Canvas
//HMS
TCanvas *hms_REF_Canv;                      //canvas to save reference time histograms
TCanvas *hhodoCanv[hod_PLANES][SIDES];
TCanvas *hcaloCanv[cal_PLANES][SIDES];
TCanvas *hdcCanv;
TCanvas *hCer_Canv;

//SHMS
TCanvas *shms_REF_Canv;                      //canvas to save reference time histograms
TCanvas *phodoCanv[hod_PLANES][SIDES];
TCanvas *pcalCanv_alt[16];
TCanvas *pPrshCanv[SIDES];
TCanvas *pdcCanv;
TCanvas *pngCer_Canv;
TCanvas *phgCer_Canv;

//TRG
TCanvas *pTRG_Canv;



// variable for setting fit limits
int binmax;
double upperlim;

double stdev;
double amplitude;
double center;
double xmin_fit; 
double xmax_fit;

TF1 * gaus_fit = 0;
  
//=========================================
//Define TLines TO DRAW AROUND CUT REGION
//=========================================

//----Reference Time TLines----

//HMS
TLine *hT1_Line;      //hms trigger ref. time
TLine *hDCREF_Line;  //hms DC ref. time
TLine *hFADC_Line;    //flash ADC ref. time
//SHMS
TLine *pT2_Line;      //shms trigger ref. time
TLine *pDCREF_Line;  //shms DC ref. time
TLine *pFADC_Line;    //flash ADC ref. time

  
//-----Detectors Time Window CUts Lines-----

//HMS
TLine *hhod_LineMin[hod_PLANES][SIDES][16];
TLine *hhod_LineMax[hod_PLANES][SIDES][16];

TLine *hcal_LineMin[cal_PLANES][SIDES][13];
TLine *hcal_LineMax[cal_PLANES][SIDES][13];

TLine *hdc_LineMin[dc_PLANES];
TLine *hdc_LineMax[dc_PLANES];

TLine *hCER_LineMin[2];
TLine *hCER_LineMax[2];

//SHMS
TLine *phod_LineMin[hod_PLANES][SIDES][21];
TLine *phod_LineMax[hod_PLANES][SIDES][21];

TLine *pcal_LineMin[224];
TLine *pcal_LineMax[224];

TLine *pPrsh_LineMin[2][14];
TLine *pPrsh_LineMax[2][14];

TLine *pdc_LineMin[dc_PLANES];
TLine *pdc_LineMax[dc_PLANES];

TLine *phgcer_LineMin[4];
TLine *phgcer_LineMax[4];

TLine *pngcer_LineMin[4];
TLine *pngcer_LineMax[4];

//TRG
TLine *ptrg1r1_LineMin;
TLine *ptrg1r1_LineMax;

TLine *ptrg1r2_LineMin;
TLine *ptrg1r2_LineMax;

TLine *ptrg4r1_LineMin;
TLine *ptrg4r1_LineMax;

TLine *ptrg4r2_LineMin;
TLine *ptrg4r2_LineMax;

//===========================================

  
//========================
//Define TTree Leaf Names 
//========================


                                                                                                                                  
//HMS Detector Leaf Names                                                                                                                                                                     
TString n_hhod_TdcAdcTimeDiff;                                                                                                                                             
TString n_hhod_AdcMult;                                                                                                                                                   
TString n_hcal_TdcAdcTimeDiff;                                                                                                                                               
TString n_hcal_AdcMult;                                                                                                                                                         
TString n_hcer_TdcAdcTimeDiff;                                                                                                                      
TString n_hcer_AdcMult;                                                                                                                                                                
TString n_hndata_rawTDC;                                                                                                                          
TString n_hdc_rawTDC; 
TString n_hdc_TdcMult;

//HMS Ref. Time Names
TString n_hT1_ref;
TString n_hDC_ref;
TString n_hFADC_ref;
TString n_hT1_tdcMult;
TString n_hDC_tdcMult;
TString n_hFADC_adcMult;

//SHMS Detector Leaf Names              
TString n_phod_TdcAdcTimeDiff;                                                                                                                                                         
TString n_phod_AdcMult;                                                                                                                                                              
TString n_pcal_TdcAdcTimeDiff;                                                                                                                                                        
TString n_pcal_AdcMult;                                                                                                                                                               
TString n_pPrSh_TdcAdcTimeDiff;
TString n_pPrSh_AdcMult;
TString n_phgcer_TdcAdcTimeDiff;                                                                                                                                          
TString n_phgcer_AdcMult;
TString n_pngcer_TdcAdcTimeDiff;                                                                                                                                                      
TString n_pngcer_AdcMult;
TString n_pdc_rawTDC;                                                                                                                                                               
TString n_pndata_rawTDC;                                                                                                                                                               
TString n_pdc_TdcMult;

//SHMS Ref. Time Names                                                                                                                                                                      
TString n_pT1_ref;                                                                                                                                                
TString n_pT2_ref;                                                                                                                                                
TString n_pDC_ref;                                                                                                                                                    
TString n_pFADC_ref;        
TString n_pT1_tdcMult;
TString n_pT2_tdcMult;                                                                                                                                                
TString n_pDC_tdcMult;                                                                                                                                                    
TString n_pFADC_adcMult;                                                                                                                                 

//TRG Detector  Leaf Names
TString n_ptrg1_r1;
TString n_ptrg1_r2;
TString n_ptrg4_r1;
TString n_ptrg4_r2;

//========================================
//Define Variables Associated with Leafs
//========================================

//HMS Leaf Variables
Double_t hhod_TdcAdcTimeDiff[hod_PLANES][SIDES][16];
Double_t hhod_AdcMult[hod_PLANES][SIDES][16];
Double_t hcer_TdcAdcTimeDiff[2];
Double_t hcer_AdcMult[2];
Double_t hcal_TdcAdcTimeDiff[cal_PLANES][SIDES][13];
Double_t hcal_AdcMult[cal_PLANES][SIDES][13];

//HMS Ref. Time Varables                                                                                                                                           
Double_t hT1_ref;                                                                                                    
Double_t hDC_ref[4];                                                                                      
Double_t hFADC_ref;
Double_t hT1_tdcMult;
Double_t hDC_tdcMult[4];
Double_t hFADC_adcMult;

//SHMS Leaf Variables
Double_t phod_TdcAdcTimeDiff[hod_PLANES][SIDES][21];
Double_t phod_AdcMult[hod_PLANES][SIDES][21];
Double_t pcal_TdcAdcTimeDiff[1][224];
Double_t pcal_AdcMult[1][224];
Double_t pPrSh_TdcAdcTimeDiff[1][SIDES][14]; 
Double_t pPrSh_AdcMult[1][SIDES][14]; 
Double_t phgcer_TdcAdcTimeDiff[4];
Double_t phgcer_AdcMult[4];
Double_t pngcer_TdcAdcTimeDiff[4];
Double_t pngcer_AdcMult[4];

//SHMS Ref. Time Varables                                                                                                                                           
Double_t pT1_ref;
Double_t pT2_ref;                                                                                                    
Double_t pDC_ref[10];                                                                                      
Double_t pFADC_ref;
Double_t pT1_tdcMult;
Double_t pT2_tdcMult;
Double_t pDC_tdcMult[10];
Double_t pFADC_adcMult;

//Drift Chamber rawTDC / Ndata / Multiplicity
Double_t hdc_rawTDC[dc_PLANES][1000];
Double_t pdc_rawTDC[dc_PLANES][1000];

Int_t hndata_rawTDC[dc_PLANES];
Int_t pndata_rawTDC[dc_PLANES];

//TRG Detector Leaf Variables
Double_t ptrg1_r1;
Double_t ptrg1_r2;
Double_t ptrg4_r1;
Double_t ptrg4_r2;

//=========================================================


#endif 
