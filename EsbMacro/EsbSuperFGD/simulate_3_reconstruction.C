/*

  Macro that shows how to run the WC digitization on an existing simulation
  (e.g., produced with ess_sim.C).

  Based on the example in the presentation from
  Konstantin Gertsenberger
  
  .L ess_dig_fgd.C
  ess_dig_fgd()
  
*/

void simulate_3_reconstruction(TString inFile = "fgd_dig.root", 
	      TString parFile = "params.root",
	      TString outFile = "fgd_recon.root",
              Int_t nStartEvent = 7, Int_t nEvents = 1)
{
  using namespace esbroot;

  FairRunAna *fRun= new FairRunAna();
  // Set Input Source and Output file
  FairFileSource *fFileSource = new FairFileSource(inFile);
  fRun->SetSource(fFileSource);

  fRun->SetSink(new FairRootFileSink(outFile));

  // -----  Parameter database   --------------------------------------------
  FairRuntimeDb* rtdb = fRun->GetRuntimeDb();

  FairParRootFileIo* parIo1 = new FairParRootFileIo();
  parIo1->open(parFile);
  rtdb->setFirstInput(parIo1);
  
  rtdb->setOutput(parIo1); 
  rtdb->saveOutput();

  double debugLvl = 0.0; 

  fair::Logger::SetConsoleSeverity(fair::Severity::debug);
  fair::Logger::SetConsoleColor(true);

  FairTask* recon = new reconstruction::superfgd::FgdGenFitRecon(
    "Reconstruction Task"             // name of the task
    ,"../../EsbGeometry/EsbSuperFGD/EsbConfig/fgdconfig"  //File with detector configuration
    ,"../../EsbGeometry/EsbSuperFGD/EsbConfig/graphConfig" // File containing graph algorithm info
    ,"../../geometry/media.geo"       // Media file with defined materials
    , 1                               // Verbose level
    , debugLvl                        // debug level of genfit (0 - little, 1 - debug info, 2 - detailed)
    , false                           // To visualize the tracks using genfit::Eventdisplay
    , "D");                           // Option to be passed for genfit::Eventdisplay if used

  ((reconstruction::superfgd::FgdGenFitRecon*)recon)->SetMinHits(3);
  // ((reconstruction::superfgd::FgdGenFitRecon*)recon)->SetUseTracker(reconstruction::superfgd::FgdGenFitRecon::TrackFinder::HOUGH_PATHFINDER_LINE);
  // ((reconstruction::superfgd::FgdGenFitRecon*)recon)->SetUseTracker(reconstruction::superfgd::FgdGenFitRecon::TrackFinder::HOUGH_PATHFINDER_HELIX);
  // ((reconstruction::superfgd::FgdGenFitRecon*)recon)->SetUseTracker(reconstruction::superfgd::FgdGenFitRecon::TrackFinder::HOUGH_PATHFINDER_CURL);
  // ((reconstruction::superfgd::FgdGenFitRecon*)recon)->SetUseTracker(reconstruction::superfgd::FgdGenFitRecon::TrackFinder::GRAPH);
  // ((reconstruction::superfgd::FgdGenFitRecon*)recon)->SetUseTracker(reconstruction::superfgd::FgdGenFitRecon::TrackFinder::GRAPH_HOUGH_PATHFINDER);
  ((reconstruction::superfgd::FgdGenFitRecon*)recon)->SetUseTracker(reconstruction::superfgd::FgdGenFitRecon::TrackFinder::HOUGH_GRAPH_LEAVES);
  

  ((reconstruction::superfgd::FgdGenFitRecon*)recon)->AddPdgMomLoss(11, 9. , 6.);
  ((reconstruction::superfgd::FgdGenFitRecon*)recon)->AddPdgMomLoss(13, 21. , 5.);
  ((reconstruction::superfgd::FgdGenFitRecon*)recon)->AddPdgMomLoss(211, 28. , 1.);
  ((reconstruction::superfgd::FgdGenFitRecon*)recon)->AddPdgMomLoss(2212, 148. , 63.);
  
  fRun->AddTask(recon);   
  fRun->Init(); // initializing
  fRun->Run(nStartEvent, nStartEvent + nEvents);
  fRun->CreateGeometryFile("geo_recon.root");  // for additional full geometry file
}
