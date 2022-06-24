#include <chrono>
#include <iostream>

using namespace std;
using namespace std::chrono;


void replay_cafe_scalers(Int_t RunNumber = 0, Int_t MaxEvent = 0, TString ftype="") {

  
  // Get RunNumber and MaxEvent if not provided.
  if(RunNumber == 0) {
    cout << "Enter a Run Number (-1 to exit): ";
    cin >> RunNumber;
    if( RunNumber<=0 ) return;
  }
  if(MaxEvent == 0) {
    cout << "\nNumber of Events to analyze: ";
    cin >> MaxEvent;
    if(MaxEvent == 0) {
      cerr << "...Invalid entry\n";
      exit;
    }
  }
  
  if(ftype==""){
    cout  << "\nEnter file type to use (e.g., shms50k, hms50k, sample, prod): \n " << endl;
    cin >> ftype;
    if(ftype==""){
      cerr << "...Invalid file type\n";
      exit;
    }
  }

  // Get starting timepoint
  auto start = high_resolution_clock::now();
    
  // Create file name patterns.
  //const char* RunFileNamePattern = "coin_all_%05d.dat";
  const char* RunFileNamePattern = "shms_all_%05d.dat";

  vector<TString> pathList;
  pathList.push_back(".");
  pathList.push_back("./raw");
  pathList.push_back("./cache");
  pathList.push_back("./raw.copiedtotape");

  //const char* RunFileNamePattern = "raw/coin_all_%05d.dat";

  // Create dir. to store monitoring histos
  TString cmd = Form("mkdir -p ROOTfiles/%s", ftype.Data());
  gSystem->Exec(cmd); // create study type dir. if it doesn't exist

  cmd = Form("mkdir -p HISTOGRAMS/%s", ftype.Data());
  gSystem->Exec(cmd); // create study type dir. if it doesn't exist

  cmd = Form("mkdir -p HISTOGRAMS/%s/PDF", ftype.Data());
  gSystem->Exec(cmd); // create study type dir. if it doesn't exist

  cmd = Form("mkdir -p HISTOGRAMS/%s/ROOT", ftype.Data());
  gSystem->Exec(cmd); // create study type dir. if it doesn't exist
  
  const char* ROOTFileNamePattern = "ROOTfiles/%s/cafe_replay_%s_%d_%d.root";

  
  // Load global parameters
  gHcParms->Define("gen_run_number", "Run Number", RunNumber);
  gHcParms->AddString("g_ctp_database_filename", "DBASE/COIN/standard.database");
  gHcParms->Load(gHcParms->GetString("g_ctp_database_filename"), RunNumber);  // load the standard.database
  gHcParms->Load(gHcParms->GetString("g_ctp_parm_filename"));                 // load the general param
  gHcParms->Load(gHcParms->GetString("g_ctp_calib_filename"));                // load the detector calib param
  gHcParms->Load(gHcParms->GetString("g_ctp_cuts_filename"));                 // load the detector cuts param 
  gHcParms->Load(gHcParms->GetString("g_ctp_kinematics_filename"), RunNumber); // load the standard.kinematics file

  //Load params for coin. trigger configuration (up-to-date)
  gHcParms->Load("PARAM/TRIG/tcoin.param"); 

  // Load the Hall C detector map
  gHcDetectorMap = new THcDetectorMap();
  gHcDetectorMap->Load("MAPS/COIN/DETEC/coin.map");

  //use spring18 config det. map  / params
  //All runs before coin 4361 did NOT have hDCREF_5 added in ROC3
  if(RunNumber<=4361){
    
    //Load params for coin. trigger configuration (spring 2018, for testing using 2018 data) 
    //gHcParms->Load("PARAM/TRIG/archive/spring18/tcoin_spring18.param");
    
    // Load 2018 map
    gHcDetectorMap->Load("MAPS/COIN/DETEC/coin_comm18.map");
  }
  
  //=:=:=:=
  // SHMS 
  //=:=:=:=
     // Dec data
  gHaApps->Add(new Podd::DecData("D","Decoder raw data"));

  // Add event handler for scaler events
  THcScalerEvtHandler* pscaler = new THcScalerEvtHandler("P", "Hall C scaler event type 1");
  pscaler->AddEvtType(1);
  pscaler->AddEvtType(4);
  pscaler->AddEvtType(5);
  pscaler->AddEvtType(6);
  pscaler->AddEvtType(7);
  pscaler->AddEvtType(129);
  pscaler->SetDelayedType(129);
  pscaler->SetUseFirstEvent(kTRUE);
  gHaEvtHandlers->Add(pscaler);

  /*
  //Add SHMS event handler for helicity scalers
  THcHelicityScaler *phelscaler = new THcHelicityScaler("P", "Hall C helicity scaler");
  //phelscaler->SetDebugFile("PHelScaler.txt");
  phelscaler->SetROC(8);   
  phelscaler->SetUseFirstEvent(kTRUE); 
  gHaEvtHandlers->Add(phelscaler); 
  */
  
  //=:=:=
  // HMS 
  //=:=:=


  // Add event handler for scaler events
  THcScalerEvtHandler *hscaler = new THcScalerEvtHandler("H", "Hall C scaler event type 4");  
  hscaler->AddEvtType(2);
  hscaler->AddEvtType(4);
  hscaler->AddEvtType(5);
  hscaler->AddEvtType(6);
  hscaler->AddEvtType(7);
  hscaler->AddEvtType(129);
  hscaler->SetDelayedType(129);
  hscaler->SetUseFirstEvent(kTRUE);
  gHaEvtHandlers->Add(hscaler);

  /*
  // Add HMS event handler for helicity scalers                                                                                             
  THcHelicityScaler *hhelscaler = new THcHelicityScaler("H", "Hall C helicity scaler"); 
  //hhelscaler->SetDebugFile("HHelScaler.txt");                                                                
  hhelscaler->SetROC(5);
  hhelscaler->SetUseFirstEvent(kTRUE); 
  gHaEvtHandlers->Add(hhelscaler);
  */


  
  // Add event handler for prestart event 125.
  THcConfigEvtHandler* ev125 = new THcConfigEvtHandler("HC", "Config Event type 125");
  gHaEvtHandlers->Add(ev125);
  // Add event handler for EPICS events
  THaEpicsEvtHandler* hcepics = new THaEpicsEvtHandler("epics", "HC EPICS event type 180");
  gHaEvtHandlers->Add(hcepics);
 
  // Set up the analyzer - we use the standard one,
  // but this could be an experiment-specific one as well.
  // The Analyzer controls the reading of the data, executes
  // tests/cuts, loops over Acpparatus's and PhysicsModules,
  // and executes the output routines.
  THcAnalyzer* analyzer = new THcAnalyzer;

  // A simple event class to be output to the resulting tree.
  // Creating your own descendant of THaEvent is one way of
  // defining and controlling the output.
  THaEvent* event = new THaEvent;

  // Define the run(s) that we want to analyze.
  // We just set up one, but this could be many.
  THcRun* run = new THcRun( pathList, Form(RunFileNamePattern, RunNumber) );

  // Set to read in Hall C run database parameters
  run->SetRunParamClass("THcRunParameters");
  
  // Eventually need to learn to skip over, or properly analyze the pedestal events
  run->SetEventRange(1, MaxEvent); // Physics Event number, does not include scaler or control events.
  run->SetNscan(1);
  run->SetDataRequired(0x7);
  run->Print();

  // Define the analysis parameters
  TString ROOTFileName = Form(ROOTFileNamePattern, ftype.Data(), ftype.Data(), RunNumber, MaxEvent);
  analyzer->SetCountMode(2);  // 0 = counter is # of physics triggers
                              // 1 = counter is # of all decode reads
                              // 2 = counter is event number

  analyzer->SetEvent(event);

  // Set EPICS event type
  analyzer->SetEpicsEvtType(180);

  // Define crate map
  analyzer->SetCrateMapFileName("MAPS/db_cratemap.dat");

  // Define output ROOT file
  analyzer->SetOutFile(ROOTFileName.Data());

  // Define DEF-file+
  TString DefTreeFile="DEF-files/cafe_scalers.def";
  analyzer->SetOdefFile(DefTreeFile);

  // Define cuts file
  TString DefCutTreeFile;

  if(ftype.Data()=="shms50k"){
    DefCutTreeFile="DEF-files/CUTS/cafe_cuts_shms.def";
  }
  else if(ftype.Data()=="hms50k"){
    DefCutTreeFile="DEF-files/CUTS/cafe_cuts_hms.def";
  }
  else{
    DefCutTreeFile="DEF-files/CUTS/cafe_scalers_cuts.def";
  }
  
  analyzer->SetCutFile(DefCutTreeFile);  // optional

  
  // File to record accounting information for cuts
  cmd = Form("mkdir -p REPORT_OUTPUT/%s", ftype.Data());  
  gSystem->Exec(cmd); // create study type dir. if it doesn't exist
  analyzer->SetSummaryFile(Form("REPORT_OUTPUT/%s/summary_cafe_%s_%d_%d.report", ftype.Data(), ftype.Data(), RunNumber, MaxEvent));  // optional

  // Start the actual analysis.
  analyzer->Process(run);

  // Create report file from template
  // C.Y. (for now we just have 1 template file, but this can be expanded to: cafe_ftype.template (e.g. cafe_heep.template, etc)
  TString REPORT_FileName=Form("REPORT_OUTPUT/%s/cafe_%s_%d_%d.report", ftype.Data(), ftype.Data(), RunNumber, MaxEvent);
  TString TEMPLATE_FileName="TEMPLATES/cafe_scalers.template";

  //testing (to use cafe template file to mimic what the shift crew would see during CaFe data-taking, commnent the 'RunNumber<3400' below)
  //if(RunNumber<=3400){
  //  REPORT_FileName=Form("REPORT_OUTPUT/%s/deut_%s_spring18.report", ftype.Data(), ftype.Data() );
  //  TEMPLATE_FileName="TEMPLATES/deut_spring18.template";
  //}
  
  analyzer->PrintReport( TEMPLATE_FileName, REPORT_FileName );  // optional

  // Write timestamp to REPORT_FILE to keep track of run start_time
  ofstream file_object;
  file_object.open(REPORT_FileName.Data(), std::ios_base::app); 
  file_object << "" << endl;
  file_object << "start_of_run: " << (run->GetDate()).AsSQLString();
  file_object << "" << endl;
  file_object.close();
  
  // Get ending timepoint
  auto stop = high_resolution_clock::now();
  
  // Get duration. Substart timepoints to
  // get duration. To cast it to proper unit
  // use duration cast method
  auto duration = duration_cast<seconds>(stop - start);
  
  cout << "Time taken by replay_cafe: "
       << duration.count() << " sec." << endl;
  
}
