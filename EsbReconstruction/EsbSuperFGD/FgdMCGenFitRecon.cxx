#include "EsbReconstruction/EsbSuperFGD/FgdMCGenFitRecon.h"
#include "EsbReconstruction/EsbSuperFGD/FgdReconTemplate.h"
#include "EsbReconstruction/EsbSuperFGD/FgdCalorimetric.h"
#include "EsbData/EsbSuperFGD/FgdDetectorPoint.h"

// FairRoot headers
#include "FairGeoBuilder.h"
#include "FairGeoInterface.h"
#include "FairGeoLoader.h"
#include "FairGeoMedia.h"
#include "FairLogger.h"
#include <FairRootManager.h>
#include "FairVolume.h"


// Root headers
#include <TClonesArray.h>
#include <TEveManager.h>
#include <TGeoElement.h>
#include <TGeoManager.h>

// Genie headers
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"

// Genfit headers
#include "AbsBField.h"
#include "AbsMeasurement.h"
#include "ConstField.h"
#include <Exception.h>
#include <EventDisplay.h>
#include <FieldManager.h>
#include "FitStatus.h"
#include <KalmanFitterRefTrack.h>
#include <KalmanFitter.h>
#include "MaterialEffects.h"
#include "MeasuredStateOnPlane.h"
#include <PlanarMeasurement.h>
#include <RKTrackRep.h>
#include "SpacepointMeasurement.h"
#include <StateOnPlane.h>
#include "TDatabasePDG.h"
#include <TGeoMaterialInterface.h>
#include <Track.h>
#include <TrackCand.h>
#include <TrackPoint.h>
#include <TRandom.h>
#include <TVector3.h>


// PathFinder headers
#include "FinderParameter.h"
#include "HoughTrafoTrackFinder.h"
#include "TrackParameterFull.h"


// STL headers
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <memory>
#include <math.h>
#include <bits/stdc++.h>

namespace esbroot {
namespace reconstruction {
namespace superfgd {

// -----   Default constructor   -------------------------------------------
FgdMCGenFitRecon::FgdMCGenFitRecon() :
  FairTask(), fsuperFgdVol(nullptr)
  , fgdConstructor("")
  , fHitArray(nullptr)
  , isDefinedMaterials(false)
  , fDebuglvl_genfit(0)
  , fmediaFile("")
  , fTracksArray(nullptr)
  , fdisplay(nullptr)
  , isGenFitVisualization(false)
  , fGenFitVisOption("")
  , fminGenFitInterations(2)
  , fmaxGenFitIterations(4)
  , fminHits(25)
  , fMCeventData("")
  , foutputRootFile("")
  , ferrOutPath("")
  , fmomErrWidth(1)
  , fMCeventNum(0)
{ 
}
// -------------------------------------------------------------------------

// -----   Constructor   -------------------------------------------
FgdMCGenFitRecon::FgdMCGenFitRecon(const char* name
                          , const char* geoConfigFile
                          , const char* mediaFile
                          , const char* eventData
                          , Int_t verbose
                          , double debugLlv
                          , bool visualize
                          , std::string visOption) :
  FairTask(name, verbose)
  , fsuperFgdVol(nullptr)
  , fgdConstructor(geoConfigFile)
  , fHitArray(nullptr)
  , isDefinedMaterials(false)
  , fDebuglvl_genfit(debugLlv)
  , fmediaFile(mediaFile)
  , fTracksArray(nullptr)
  , fdisplay(nullptr)
  , isGenFitVisualization(visualize)
  , fGenFitVisOption(visOption)
  , fminGenFitInterations(2)
  , fmaxGenFitIterations(4)
  , fminHits(25)
  , fMCeventData(eventData)
  , foutputRootFile("")
  , ferrOutPath("")
  , fmomErrWidth(1)
  , fMCeventNum(0)
{ 
  fParams.LoadPartParams(geoConfigFile);
}
// -------------------------------------------------------------------------



// -----   Destructor   ----------------------------------------------------
FgdMCGenFitRecon::~FgdMCGenFitRecon() 
{
  if(fHitArray) {
    fHitArray->Delete();
    delete fHitArray;
  }

  if(fTracksArray) {
    fTracksArray->Delete();
    delete fTracksArray;
  }
}
// -------------------------------------------------------------------------



// -----   Public method Init   --------------------------------------------
InitStatus FgdMCGenFitRecon::Init() 
{   
  // Create the real Fgd geometry
  DefineMaterials();
  fsuperFgdVol = fgdConstructor.Construct();
  gGeoManager->SetTopVolume(fsuperFgdVol); 

  // Get dimentions from geometry file
  flunit = fParams.GetLenghtUnit(); // [cm]

  f_step_X  = fParams.ParamAsDouble(esbroot::geometry::superfgd::DP::length_X) * flunit;
  f_step_Y  = fParams.ParamAsDouble(esbroot::geometry::superfgd::DP::length_Y) * flunit;
  f_step_Z  = fParams.ParamAsDouble(esbroot::geometry::superfgd::DP::length_Z) * flunit;

  f_bin_X = fParams.ParamAsInt(esbroot::geometry::superfgd::DP::number_cubes_X);
  f_bin_Y = fParams.ParamAsInt(esbroot::geometry::superfgd::DP::number_cubes_Y);
  f_bin_Z = fParams.ParamAsInt(esbroot::geometry::superfgd::DP::number_cubes_Z);

  f_total_X = f_step_X * f_bin_X;
  f_total_Y = f_step_Y * f_bin_Y;
  f_total_Z = f_step_Z * f_bin_Z;

  //init geometry and mag. field
  TVector3 magField = fgdConstructor.GetMagneticField(); // values are in kGauss
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(magField.X(),magField.Y(), magField.Z())); 
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
  genfit::MaterialEffects::getInstance()->setEnergyLossBetheBloch(true);
  genfit::MaterialEffects::getInstance()->setEnergyLossBrems(true);
  genfit::MaterialEffects::getInstance()->setDebugLvl(fDebuglvl_genfit);

  std::ifstream eventFileStream;
  try
  {        
      eventFileStream.open(fMCeventData.c_str(), std::ios::in);

      if(eventFileStream.is_open())
      {
          std::string line;
          while(std::getline(eventFileStream,line))
          {
            fMCeventRecords.emplace_back(FgdTMVAEventRecord(line));
          }
      }
  }
  catch(const std::exception& e)
  {
      LOG(fatal) << e.what();
  }

  if(eventFileStream.is_open())
  {
      eventFileStream.close();
  } 

  // Get RootManager
  FairRootManager* manager = FairRootManager::Instance();
  if ( !manager ) {
    LOG(fatal) << "-E- FgdGenFitRecon::Init: " << "FairRootManager not instantised!";
    return kFATAL;
  }

  fHitArray = (TClonesArray*) manager->GetObject(geometry::superfgd::DP::FGD_HIT.c_str());
  if (!fHitArray) 
  {
      LOG(fatal) << "Exec", "No fgd hits array";
      return kFATAL;
  }

  OutputFileInit(manager);

  if(isGenFitVisualization)
  {
    fdisplay = genfit::EventDisplay::getInstance();
  }
  
  if(fdisplay!=nullptr && !fGenFitVisOption.empty())
  {
    fdisplay->setOptions(fGenFitVisOption);
  }

  return kSUCCESS;
}

void FgdMCGenFitRecon::OutputFileInit(FairRootManager* manager)
{
  // Create and register output array
  fTracksArray = new TClonesArray(genfit::Track::Class(), 1000);
  manager->Register(geometry::superfgd::DP::FGD_FIT_TRACK.c_str()
                    , geometry::superfgd::DP::FGD_BRANCH_FIT.c_str()
                    , fTracksArray, kTRUE);
}


// -------------------------------------------------------------------------



// -----   Public methods   --------------------------------------------
void FgdMCGenFitRecon::FinishEvent()
{
  if(isGenFitVisualization && !fgenTracks.empty())
  {
    fdisplay->addEvent(fgenTracks);
    fgenTracks.clear();
  }
  fTracksArray->Delete();
}

void FgdMCGenFitRecon::FinishTask()
{
  if(isGenFitVisualization)
  {
    fdisplay->open();
  }

  static const TVector3 z_axis(0,0,1);

  if(foutputRootFile.empty()) return;

  TFile * outFile = new TFile(foutputRootFile.c_str(), "RECREATE", "Fitted TVMA data from Fgd Detector");
  outFile->SetCompressionLevel(9);

  try
  {
    TTree * trainTree = new TTree("trainTree"
                                ,esbroot::geometry::superfgd::DP::FGD_TMVA_DATA_ROOT_FILE.c_str());

    Float_t lnuEnergy = 0.;
    Float_t muon_mom = 0;
    Float_t muon_Angle = 0;
    Float_t totPh = 0.;
    Float_t totCubes = 0;

    trainTree->Branch("totalCubes", &totCubes);
    trainTree->Branch("totalPhotons", &totPh);
    trainTree->Branch("nuEnergy", &lnuEnergy);
    trainTree->Branch("muon_mom", &muon_mom);
    trainTree->Branch("muon_angle", &muon_Angle);

    //========================================================
    TTree * fittedMomTree = new TTree("fittedMomTree"
                                ,esbroot::geometry::superfgd::DP::FGD_TMVA_DATA_ROOT_FILE.c_str());


    Float_t fit_lnuEnergy = 0.;
    Float_t fit_muon_mom = 0;
    Float_t fit_muon_Angle = 0;
    Float_t fit_totPh = 0.;
    Float_t fit_totCubes = 0;

    fittedMomTree->Branch("totalCubes", &fit_totCubes);
    fittedMomTree->Branch("totalPhotons", &fit_totPh);
    fittedMomTree->Branch("nuEnergy", &fit_lnuEnergy);
    fittedMomTree->Branch("muon_mom", &fit_muon_mom);
    fittedMomTree->Branch("muon_angle", &fit_muon_Angle);

    //========================================================
    TTree * calorimetricMomTree = new TTree("CalMomTree"
                                ,esbroot::geometry::superfgd::DP::FGD_TMVA_DATA_ROOT_FILE.c_str());


    Float_t cal_lnuEnergy = 0.;
    Float_t cal_muon_mom = 0;
    Float_t cal_muon_Angle = 0;
    Float_t cal_totPh = 0.;
    Float_t cal_totCubes = 0;

    calorimetricMomTree->Branch("totalCubes", &cal_totCubes);
    calorimetricMomTree->Branch("totalPhotons", &cal_totPh);
    calorimetricMomTree->Branch("nuEnergy", &cal_lnuEnergy);
    calorimetricMomTree->Branch("muon_mom", &cal_muon_mom);
    calorimetricMomTree->Branch("muon_angle", &cal_muon_Angle);

    const Int_t evInd = fMCeventRecords.size();
    
    Int_t limit = evInd < fMCeventNum ? evInd : fMCeventNum;

    FgdTMVAEventRecord* dataEvent = nullptr;
    for(size_t ind = 0 ; ind < limit; ind++)
    {
        LOG(debug2) << "Writing data for event " << ind;
        dataEvent = &fMCeventRecords[ind];

        bool isQuasiCC = dataEvent->IsWeakCC() && dataEvent->IsQuasiElastic();
        if(!isQuasiCC)
        {
            continue;
        }

        lnuEnergy = dataEvent->GetNuE();
        // const std::vector<std::pair<Int_t, TVector3>>& particles = dataEvent->GetPrimaryParticles();
        // for(size_t p = 0; p < particles.size(); ++p)
        // {
        //   std::pair<Int_t, TVector3> pp = particles[p];
        //   if(genie::pdg::IsMuon(pp.first) || genie::pdg::IsAntiMuon(pp.first))
        //   {
        //     muon_mom = pp.second.Mag();
        //     break;
        //   }
        // }

        muon_mom = dataEvent->GetMuonMom().Mag();
        muon_Angle = z_axis.Angle(dataEvent->GetMuonMom());
        totPh = dataEvent->GetTotalPhotons().X() + dataEvent->GetTotalPhotons().Y() + dataEvent->GetTotalPhotons().Z();
        totCubes = dataEvent->GetTotalCubes();

        fit_lnuEnergy = lnuEnergy;
        fit_totPh = totPh;
        fit_totCubes = totCubes;
        fit_muon_mom = dataEvent->GetMuonFitMom().Mag();
        fit_muon_Angle = z_axis.Angle(dataEvent->GetMuonFitMom());

        cal_lnuEnergy = lnuEnergy;
        cal_totPh = totPh;
        cal_totCubes = totCubes;
        cal_muon_mom = dataEvent->GetMuonCalorimetricMom().Mag();
        cal_muon_Angle = z_axis.Angle(dataEvent->GetMuonCalorimetricMom());

        trainTree->Fill();
        fittedMomTree->Fill();
        calorimetricMomTree->Fill();
    }

    //========================================================
    TTree * electronFitErrTree = new TTree("electronFitErrTree"
                                ,esbroot::geometry::superfgd::DP::FGD_TMVA_DATA_ROOT_FILE.c_str());

    Float_t electronFitErr = 0.;
    electronFitErrTree->Branch("electronFitErr", &electronFitErr);
    for(size_t i = 0; i < felectronFitErr.size(); ++i)
    {
      electronFitErr = felectronFitErr[i];
      electronFitErrTree->Fill();
    }
    //========================================================

    //========================================================
    TTree * muonFitErrTree = new TTree("muonFitErrTree"
                                ,esbroot::geometry::superfgd::DP::FGD_TMVA_DATA_ROOT_FILE.c_str());

    Float_t muonFitErr = 0;
    muonFitErrTree->Branch("muonFitErr", &muonFitErr);
    for(size_t i = 0; i < fmuonFitErr.size(); ++i)
    {
      muonFitErr = fmuonFitErr[i];
      muonFitErrTree->Fill();
    }
    //========================================================

    //========================================================
    TTree * protonFitErrTree = new TTree("protonFitErrTree"
                                ,esbroot::geometry::superfgd::DP::FGD_TMVA_DATA_ROOT_FILE.c_str());

    Float_t protonFitErr = 0;;
    protonFitErrTree->Branch("protonFitErr", &protonFitErr);
    for(size_t i = 0; i < fprotonFitErr.size(); ++i)
    {
      protonFitErr = fprotonFitErr[i];
      protonFitErrTree->Fill();
    }
    //========================================================

    //========================================================
    TTree * pionFitErrTree = new TTree("pionFitErrTree"
                                ,esbroot::geometry::superfgd::DP::FGD_TMVA_DATA_ROOT_FILE.c_str());

    Float_t pionFitErr = 0;;
    pionFitErrTree->Branch("pionFitErr", &pionFitErr);
    for(size_t i = 0; i < fpionFitErr.size(); ++i)
    {
      pionFitErr = fpionFitErr[i];
      pionFitErrTree->Fill();
    }
    //========================================================
  
    outFile->WriteTObject(trainTree);  
    outFile->WriteTObject(fittedMomTree);  
    outFile->WriteTObject(calorimetricMomTree);     

    outFile->WriteTObject(electronFitErrTree); 
    outFile->WriteTObject(muonFitErrTree); 
    outFile->WriteTObject(protonFitErrTree); 
    outFile->WriteTObject(pionFitErrTree);                 
  }
  catch(...)
  {

  }
  
  outFile->Close();
  
  delete outFile;


  //========================================================
  //==========  Write Err data
  //========================================================
  if(!ferrOutPath.empty())
  {
    writeErrFile(std::string("_muon"),  Axis::X,      fmuonFitData);
    writeErrFile(std::string("_muon"),  Axis::Y,      fmuonFitData);
    writeErrFile(std::string("_muon"),  Axis::Z,      fmuonFitData);

    writeErrFile(std::string("_proton"), Axis::X,     fprotonFitData);
    writeErrFile(std::string("_proton"), Axis::Y,     fprotonFitData);
    writeErrFile(std::string("_proton"), Axis::Z,     fprotonFitData);

    writeErrFile(std::string("_electron"),  Axis::X,  felectronFitData);
    writeErrFile(std::string("_electron"),  Axis::Y,  felectronFitData);
    writeErrFile(std::string("_electron"),  Axis::Z,  felectronFitData);

    writeErrFile(std::string("_pion"),  Axis::X,      fpionFitData);
    writeErrFile(std::string("_pion"),  Axis::Y,      fpionFitData);
    writeErrFile(std::string("_pion"),  Axis::Z,      fpionFitData);
  }
}

void FgdMCGenFitRecon::writeErrFile(const std::string& fileEnding, Axis axis, std::vector<FitData>& dataVec)
{
  std::stringstream ss;
  ss << ferrOutPath << fileEnding;
  switch(axis)
  {
    case Axis::X: ss << "_X";
          break;
    case Axis::Y: ss << "_Y";
          break;
    case Axis::Z: ss << "_Z";
          break;
    default: 
          break;
  }
  ss<< ".txt";
  std::string path = ss.str();

  std::map<int, std::vector<float>> projErr;
  for(size_t i = 0; i < dataVec.size(); ++i)
  {
    FitData& data = dataVec[i];
    int proj(0);
    float err(0.);
    switch(axis)
    {
      case Axis::X: proj = data.trackProj.X();
                    err = (data.mcMom.X() - data.fitMom.X());
                    break;
      case Axis::Y: proj = data.trackProj.Y();
                    err = (data.mcMom.Y() - data.fitMom.Y());
                    break;
      case Axis::Z: proj = data.trackProj.Z(); 
                    err = (data.mcMom.Z() - data.fitMom.Z());;
                    break;
      default: 
                    break;
    }

    bool hasProj = (projErr.find(proj) != projErr.end());
    if(hasProj)
    {
      projErr[proj].insert(projErr[proj].end(), err);
    }
    else
    {
      std::vector<float> errs;
      errs.insert(errs.end(), err);
      projErr[proj] = std::move(errs);
    }
  }

  std::ofstream outfile(path, std::ios::trunc);
  if(outfile.is_open())
  {
    auto it = projErr.begin();
    while(it!= projErr.end())
    {
      int keyProj = it->first;
      float stdDev = calcStdDev(it->second);
      ++it;

      outfile << keyProj << " " << stdDev << endl;
    }
  }
  outfile.close();
}

float  FgdMCGenFitRecon::calcStdDev(std::vector<float>& vals)
{
  float sum(0.);
  for(size_t i = 0; i < vals.size(); ++i)
  {
    sum+= vals[i];
  }
  float mean = sum/vals.size();

  float sumMean(0.);
  for(size_t i = 0; i < vals.size(); ++i)
  {
    float&& temp = (vals[i] - mean);
    sumMean+= (temp*temp);
  }

  float avgsum = sumMean/vals.size();
  return std::sqrt(avgsum);
}

void FgdMCGenFitRecon::Exec(Option_t* opt) 
{  
  try
  {
    std::vector<ReconHit> allhits;
    std::vector<std::vector<ReconHit>> foundTracks;

    bool rc = GetHits(allhits);

    if(rc)
    { 
      LOG(debug) <<" Hits to retrieve tracks from " << allhits.size();
      SplitTrack(allhits, foundTracks);
    }

    if(rc)
    {
      LOG(debug) <<" Found tracks to fit " << foundTracks.size();
      FitTracks(foundTracks);
    }
    else
    {
      LOG(error) << " Could not find clean hits or tracks! ";
    }
  }
  catch(genfit::Exception& e)
  {
      LOG(fatal) << "Exception, when tryng to fit track";
      LOG(fatal) << e.what();
  }
  ++fMCeventNum;
}
// -------------------------------------------------------------------------


// -----   Protected methods   --------------------------------------------

Bool_t FgdMCGenFitRecon::GetHits(std::vector<ReconHit>& allHits)
{
  Double_t errPhotoLimit = fParams.ParamAsDouble(esbroot::geometry::superfgd::DP::FGD_ERR_PHOTO_LIMIT);

  LOG(debug) << "fHitArray->GetEntries() " << fHitArray->GetEntries();

  std::map<Long_t, Bool_t> visited;

  Int_t numVis(0);
 
  for(Int_t i =0; i <  fHitArray->GetEntries() ; i++)
  {
    data::superfgd::FgdHit* hit = (data::superfgd::FgdHit*)fHitArray->At(i);
    TVector3  photoE = hit->GetPhotoE();    
    TVector3  mppcLoc = hit->GetMppcLoc();  

    Int_t&& x = mppcLoc.X();
    Int_t&& y = mppcLoc.Y();
    Int_t&& z = mppcLoc.Z();

    // LOG(debug2) << "TrackId " << hit->GetTrackId();
    // LOG(debug2) << "GetPgd " << hit->GetPgd();
    // LOG(debug2) << "GetTrackLengthOrigin " << hit->GetTrackLengthOrigin();
    // LOG(debug2) << "GetTrackLenght " << hit->GetTrackLenght();
    // LOG(debug2) << "Photons " <<"(errLimit = " << errPhotoLimit << ") " << photoE.X() << " " << photoE.Y()<< " " << photoE.Z();
    // LOG(debug2) << "Pos "  << x << " " << y << " " << z;
    
    // Int_t ind = ArrInd(x,y,z);
    // LOG(debug2) << "ArrInd " << ind;
    // if(visited[ind])
    // {
    //   // TODO: fix, causes "free(): invalid next size (fast)""
    //   // If already exists, add the photons
    //   // ReconHit toFind;
    //   // toFind.fmppcLoc = mppcLoc;
    //   // std::vector<ReconHit>::iterator recHit = find(allHits.begin(), allHits.end(), toFind);
    //   // ReconHit& foundHit = *recHit;
    //   // foundHit.fphotons = foundHit.fphotons + photoE;
    //   LOG(debug) << "isvisited ";
    //   continue;
    // }
    // visited[ind] = true;

    Double_t totalPhotons = photoE.X() + photoE.Y() + photoE.Z();
    if(totalPhotons >= errPhotoLimit)
    {
      TVectorD hitPos(3);
      hitPos(0) = -f_total_X/2 + f_step_X*mppcLoc.X()  +f_step_X/2;
      hitPos(1) = -f_total_Y/2 + f_step_Y*mppcLoc.Y()  +f_step_Y/2;
      hitPos(2) = -f_total_Z/2 + f_step_Z*mppcLoc.Z()  +f_step_Z/2;

      allHits.emplace_back(ReconHit(
                                mppcLoc
                              , TVector3(hitPos(0),hitPos(1),hitPos(2))
                              , hit->GetDpos()
                              , photoE
                              , hit->GetTime()
                              , hit->GetMomentum()
                              , hit->GetExitMomentum()
                              , hit->GetTrackLenght()
                              , hit->GetTrackLengthOrigin()
                              , hit->GetPgd()
                              , hit->GetTrackId()
                              , hit->GetEdep()
                              , hit->GetPhotoDist1()
                              , hit->GetMppcDist1()
                              , hit->GetPhotoDist2()
                              , hit->GetMppcDist2()
                              , hit->GetPe()
                            )); 
    }
  }

  LOG(debug) << "allHits.size() " << allHits.size();

  return (allHits.size() > 0);
}


void FgdMCGenFitRecon::SplitTrack(std::vector<ReconHit>& allHits, std::vector<std::vector<ReconHit>>& foundTracks)
{
  std::map<Int_t, std::vector<ReconHit>> tracks;
  for(size_t i = 0; i < allHits.size(); ++i)
  {
      ReconHit& rh = allHits[i];

      if(tracks.find(rh.ftrackId)!=tracks.end())
      {
          tracks[rh.ftrackId].push_back(rh);
      }
      else
      {
          tracks[rh.ftrackId] = std::vector<ReconHit> {rh};
          LOG(debug2) << rh.ftrackId << " " << rh.fpdg << "[pdg]";
      }
  }

  for(auto iter = tracks.begin(); iter!=tracks.end(); ++iter)
  {
      LOG(debug2) << iter->first << " track size " << (iter->second).size();
      foundTracks.emplace_back(iter->second);
  }

  LOG(info) << "Found tracks " << foundTracks.size();
}


void FgdMCGenFitRecon::FitTracks(std::vector<std::vector<ReconHit>>& foundTracks)
{
    // init fitter
    std::shared_ptr<genfit::AbsKalmanFitter> fitter = make_shared<genfit::KalmanFitterRefTrack>();
    //std::shared_ptr<genfit::AbsKalmanFitter> fitter = make_shared<genfit::KalmanFitter>();
    fitter->setMinIterations(fminGenFitInterations);
    fitter->setMaxIterations(fmaxGenFitIterations);
    //fitter->setDeltaPval(1e-7);
    fitter->setDebugLvl(fDebuglvl_genfit);

    bool primaryMuonFound(false);

    int detId(1); // Detector id, it is the same, we only have one detector

    for(size_t i = 0; i <  foundTracks.size() ; ++i)
    {
      std::vector<ReconHit>& hitsOnTrack = foundTracks[i];
      // Sort by time, the 1st hit in time is the start of the track
      std::sort(hitsOnTrack.begin(), hitsOnTrack.end(), [](ReconHit& bh1, ReconHit& bh2){return bh1.ftime<bh2.ftime;});
    }

    FgdTMVAEventRecord& tvmaEventRecord = fMCeventRecords[fMCeventNum];

    Int_t sumTotalCubes = 0;
    TVector3 sumTotalPhoto(0,0,0);

    // Extract TVMA photon and cube data
    for(size_t i = 0; i <  foundTracks.size() ; ++i)
    {
        std::vector<ReconHit>& track = foundTracks[i];
        if(isParticleNeutral(track[0].fpdg))
        {
          continue;
        }
        for(size_t j = 0; j <  track.size() ; ++j)
        {
            ReconHit& hit = track[j];
            sumTotalPhoto +=hit.fphotons; 
            ++sumTotalCubes;
        } 
    }

    tvmaEventRecord.SetTotalPhotons(sumTotalPhoto);
    tvmaEventRecord.SetTotalCubes(sumTotalCubes);
    // =========================================

    for(size_t i = 0; i <  foundTracks.size() ; ++i)
    {
      std::vector<ReconHit>& hitsOnTrack = foundTracks[i];
      bool isMuontrack = (hitsOnTrack[0].fpdg == genie::kPdgMuon || hitsOnTrack[0].fpdg == genie::kPdgAntiMuon);

      Int_t maxPhotonsPerCube = 1;
      for(Int_t ph = 0; ph < hitsOnTrack.size(); ++ph)
      {
        Int_t sumPhotons = hitsOnTrack[ph].fphotons.X() + hitsOnTrack[ph].fphotons.Y() + hitsOnTrack[ph].fphotons.Z();
        if(maxPhotonsPerCube < sumPhotons) {maxPhotonsPerCube = sumPhotons;}
      }

      // Set lower limit on track size
      if(hitsOnTrack.size()<fminHits)
      {
        LOG(debug) << "Track " << i << " below limit, continue with next track (" << hitsOnTrack.size() << " < " << fminHits << ")";
        continue;
      }

      if(!isAllowed(hitsOnTrack[0].fpdg) || isParticleNeutral(hitsOnTrack[0].fpdg))
      {
        LOG(debug) << " Particle is not allowed for fitting " << hitsOnTrack[0].fpdg;
        continue;
      }
      LOG(debug) << " Fit particle track " << hitsOnTrack[0].ftrackId << " " 
              <<  hitsOnTrack[0].fpdg << " [Momentum = " << hitsOnTrack[0].fmom.Mag() << "]"
              << "(" << hitsOnTrack[0].fmom.X() << "," << hitsOnTrack[0].fmom.Y() << "," << hitsOnTrack[0].fmom.Z() << ")";
      
      // if(!isMuontrack)
      // {
      //   LOG(debug) << " Fitting only muon tracks pdg = " << hitsOnTrack[0].fpdg;
      //   continue;
      // }



      double sum = 0.;
      for(Int_t bh = 0; bh < hitsOnTrack.size(); ++bh)
      {
        std::sort(hitsOnTrack.begin(), hitsOnTrack.end(), [](ReconHit& bh1, ReconHit& bh2){return bh1.ftime<bh2.ftime;});
        ReconHit& hit = hitsOnTrack[bh];
        sum += hit.fEdep;
        LOG(warning) << " \tMom [ "  << hit.fpdg << " ] \t" << hit.fmom.Mag()  << "  \tdEdx " << hit.fEdep
                    << "  \tL " << hit.ftrackLength << "  \tE/L " << hit.fEdep/hit.ftrackLength
                    << "(" << hit.fmppcLoc.X() << "," << hit.fmppcLoc.Y() << "," << hit.fmppcLoc.Z() << ")";
      }

      LOG(warning) << " \tAvg dEdx " << 12.71 * sum/hitsOnTrack[hitsOnTrack.size()-1].ftrackLengthOrigin;

      //continue;

      // Sort by time, the 1st hit in time is the start of the track
      //std::sort(hitsOnTrack.begin(), hitsOnTrack.end(), [](ReconHit& bh1, ReconHit& bh2){return bh1.ftime<bh2.ftime;});
      
      const int pdg = hitsOnTrack[0].fpdg;
      TVector3 posM(hitsOnTrack[0].fHitPos);
      TVector3 momM(hitsOnTrack[0].fmom);

      TVector3 calMom = getCalorimetricMomentum(hitsOnTrack);
      // momM = calMom;
      // if(isMuontrack)
      // {
      //     // calMom = getCalorimetricMomentum(hitsOnTrack);
      //     // momM = calMom;
      //     momM = tvmaEventRecord.GetMuonMom();
      // }

      //Bool_t inlcudeBetheBloch = (pdg != genie::kPdgProton); // For proton ignore bethe bloch
      //genfit::MaterialEffects::getInstance()->setEnergyLossBetheBloch(inlcudeBetheBloch);

      // genfit::MaterialEffects::getInstance()->setEnergyLossBetheBloch(true);
      // genfit::MaterialEffects::getInstance()->setMscModel("Highland");
      // genfit::MaterialEffects::getInstance()->drawdEdx(pdg);
      // genfit::MaterialEffects::getInstance()->setNoEffects();
      // genfit::MaterialEffects::getInstance()->setNoiseBetheBloch(true);
      // genfit::MaterialEffects::getInstance()->setNoiseCoulomb(true);
      // genfit::MaterialEffects::getInstance()->setEnergyLossBrems(true);
      // genfit::MaterialEffects::getInstance()->setNoiseBrems(true);
      // genfit::MaterialEffects::getInstance()->ignoreBoundariesBetweenEqualMaterials(true);


      // if(isParticleNeutral(pdg))
      // {
      //   LOG(debug) << "Track " << i << " is of neutral particle ["<< pdg << "] continue with next track.";
      //   continue;
      // }

      // approximate covariance
      static const double resolution = 0.1;// Default in example is 0.1;   1/sqrt(12)
                                    // cm; resolution of generated measurements
      TMatrixDSym hitCov(3);
      hitCov(0,0) = resolution*resolution;
      hitCov(1,1) = resolution*resolution;
      hitCov(2,2) = resolution*resolution;

      TMatrixDSym covM(6);
      for (int ci = 0; ci < 3; ++ci)
          covM(ci,ci) = resolution*resolution;
      for (int ci = 3; ci < 6; ++ci)
          covM(ci,ci) = covM(ci,ci) = pow(  ((resolution / hitsOnTrack.size()) / sqrt(3)), 2); 

      // trackrep
      genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);

      // smeared start state
      genfit::MeasuredStateOnPlane stateSmeared(rep);
      stateSmeared.setPosMomCov(posM, momM, covM);
      

      // create track
      TVectorD seedState(6);
      TMatrixDSym seedCov(6);
      stateSmeared.get6DStateCov(seedState, seedCov);
  
      //genfit::Track* toFitTrack = new genfit::Track(rep, seedState, seedCov);
      genfit::Track* toFitTrack = new genfit::Track(rep, posM, momM);

      LOG(debug) << "******************************************* ";
      LOG(debug) << "******    Track "<< i << "  ************************";
      LOG(debug) << "******************************************* ";
      LOG(debug) << " \tPdg code " << pdg;
      LOG(debug) << " \tHits in track "<< hitsOnTrack.size();
      LOG(debug) << " \tTrack Momentum [" << hitsOnTrack[0].fmom.Mag() << "]" << "(" << hitsOnTrack[0].fmom.X() << "," << hitsOnTrack[0].fmom.Y() << "," << hitsOnTrack[0].fmom.Z() << ")";

      // if(isMuontrack)
      // {
      //     LOG(debug) << " \tCalorimetric Momentum [" << calMom.Mag() << "]" << "(" << calMom.X() << "," << calMom.Y() << "," << calMom.Z() << ")";
      // }

      //for(Int_t bh = 0; bh < hitsOnTrack.size(); bh+=3)
      //std::vector<genfit::AbsMeasurement*> measurements;
      for(Int_t bh = 0; bh < hitsOnTrack.size(); ++bh)
      {
        TVectorD hitPos(3);
        hitPos(0) = hitsOnTrack[bh].fHitPos.X();
        hitPos(1) = hitsOnTrack[bh].fHitPos.Y();
        hitPos(2) = hitsOnTrack[bh].fHitPos.Z();

        // hitPos(0) = hitsOnTrack[bh].fMCPos.X();
        // hitPos(1) = hitsOnTrack[bh].fMCPos.Y();
        // hitPos(2) = hitsOnTrack[bh].fMCPos.Z();
        // LOG(warning) << "X " << hitsOnTrack[bh].fHitPos.X() << " " << hitsOnTrack[bh].fMCPos.X();
        // LOG(warning) << "Y " << hitsOnTrack[bh].fHitPos.Y() << " " << hitsOnTrack[bh].fMCPos.Y();
        // LOG(warning) << "Z " << hitsOnTrack[bh].fHitPos.Z() << " " << hitsOnTrack[bh].fMCPos.Z();

        
        Int_t sumPhotons = hitsOnTrack[bh].fphotons.X() + hitsOnTrack[bh].fphotons.Y() + hitsOnTrack[bh].fphotons.Z();

        TMatrixDSym phhitCov(3);
        // The more photons, the more likely the particle went through the center
        // dimention of resolution is in [cm]
        double cube_diag_length = (sqrt(3)/2);
        double normCoeff = cube_diag_length*(1 - sumPhotons/maxPhotonsPerCube);     //std::sqrt(sumPhotons);
        //normCoeff = std::sqrt(normCoeff);
        normCoeff = normCoeff == 0 ? 1 : normCoeff;
        phhitCov(0,0) = resolution * resolution  * normCoeff;
        phhitCov(1,1) = resolution * resolution  * normCoeff;
        phhitCov(2,2) = resolution * resolution  * normCoeff;


        //genfit::AbsMeasurement* measurement = new genfit::SpacepointMeasurement(hitPos, hitCov, detId, bh, nullptr);
        //genfit::AbsMeasurement* measurement = new genfit::SpacepointMeasurement(hitPos, phhitCov, detId, bh, nullptr);
        //measurements.emplace_back(measurement);

        genfit::AbsMeasurement* measurement = new genfit::SpacepointMeasurement(hitPos, phhitCov, detId, bh, nullptr);
        toFitTrack->insertMeasurement(measurement, bh);
      }
      //toFitTrack->insertPoint(new genfit::TrackPoint(measurements, toFitTrack));


      // for(Int_t bh = 0; bh < hitsOnTrack.size()-1; ++bh)
      // {
      //   ReconHit& first = hitsOnTrack[bh];
      //   ReconHit& second = hitsOnTrack[bh+1];
      //   TVectorD hitPos(3);
      //   hitPos(0) = (first.fHitPos.X() + second.fHitPos.X())/2;
      //   hitPos(1) = (first.fHitPos.Y() + second.fHitPos.Y())/2;
      //   hitPos(2) = (first.fHitPos.Z() + second.fHitPos.Z())/2;;
        

      //   genfit::AbsMeasurement* measurement = new genfit::SpacepointMeasurement(hitPos, hitCov, detId, bh + hitsOnTrack.size(), nullptr);
      //   //measurements.emplace_back(measurement);
      //   toFitTrack->insertMeasurement(measurement, bh);
      // }


      try
      {
        //check
        toFitTrack->checkConsistency();

        // do the fit
        fitter->processTrack(toFitTrack, true);

        //check
        toFitTrack->checkConsistency();

        PrintFitTrack(*toFitTrack);
        genfit::FitStatus* fiStatuStatus = toFitTrack->getFitStatus();

        WriteOutput(  pdg
                    , (*toFitTrack).getFittedState().getMom()
                    , hitsOnTrack[0].fmom
                    , *toFitTrack
                    , fiStatuStatus
                    , hitsOnTrack);

        if(fiStatuStatus->isFitted() && fiStatuStatus->isFitConverged())
        {
          fgenTracks.push_back(toFitTrack);
        }

        if(isMuontrack && !primaryMuonFound )
        {
          if(fiStatuStatus->isFitted() && (fiStatuStatus->isFitConverged() || fiStatuStatus->isFitConvergedPartially()) )
          {
            TVector3 fitmom = (toFitTrack->getFittedState()).getMom();
            tvmaEventRecord.SetMuonFitMom(fitmom);
            primaryMuonFound = true;
          }

          tvmaEventRecord.SetMuonCalorimetricMom(calMom);
        }
        
        LOG(debug) <<"******************************************* ";
      }
      catch(genfit::Exception& e)
      {
          LOG(error) <<"Exception, when tryng to fit track";
          LOG(error) << e.what();
          LOG(error) << e.getExcString();

          tvmaEventRecord.SetMuonCalorimetricMom(calMom); 
      }
    }
}

void FgdMCGenFitRecon::DefineMaterials() 
{
  if(isDefinedMaterials) return; // Define materials only once

  isDefinedMaterials = true;

  new FairGeoLoader("TGeo","Geo Loader");
  FairGeoLoader *geoLoad = FairGeoLoader::Instance();
  if(geoLoad==nullptr)
  {
    LOG(error)<< "geoLoad is null";
    std::cout << "geoLoad is null" << endl;
    throw;
  }

	FairGeoInterface *geoFace = geoLoad->getGeoInterface();

  geoFace->setMediaFile(fmediaFile.c_str());
  geoFace->readMedia();

	FairGeoMedia *geoMedia = geoFace->getMedia();
	FairGeoBuilder* geoBuild = geoLoad->getGeoBuilder();

  // FairGeoMedium* brass = geoMedia->getMedium(esbroot::geometry::superfgd::materials::brass);
	// geoBuild->createMedium(brass);

  // FairGeoMedium* bronze = geoMedia->getMedium(esbroot::geometry::superfgd::materials::bronze);
	// geoBuild->createMedium(bronze);

  // FairGeoMedium* stainlessSteel = geoMedia->getMedium(esbroot::geometry::superfgd::materials::stainlessSteel);
	// geoBuild->createMedium(stainlessSteel);

  // FairGeoMedium* methane = geoMedia->getMedium(esbroot::geometry::superfgd::materials::methane);
	// geoBuild->createMedium(methane);

  // FairGeoMedium* carbonDioxide = geoMedia->getMedium(esbroot::geometry::superfgd::materials::carbonDioxide);
	// geoBuild->createMedium(carbonDioxide);

  // FairGeoMedium* carbontetraFloride = geoMedia->getMedium(esbroot::geometry::superfgd::materials::carbontetraFloride);
	// geoBuild->createMedium(carbontetraFloride);

  // FairGeoMedium* titaniumDioxide = geoMedia->getMedium(esbroot::geometry::superfgd::materials::titaniumDioxide);
	// geoBuild->createMedium(titaniumDioxide);

  // FairGeoMedium* polystyrene = geoMedia->getMedium(esbroot::geometry::superfgd::materials::polystyrene);
	// geoBuild->createMedium(polystyrene);

  FairGeoMedium* scintillator = geoMedia->getMedium(esbroot::geometry::superfgd::materials::scintillator);
  scintillator->setMediumIndex(esbroot::geometry::superfgd::materials::GetNextIndex());
	geoBuild->createMedium(scintillator);
  scintillator->Print();

  FairGeoMedium* paraterphnyl = geoMedia->getMedium(esbroot::geometry::superfgd::materials::paraterphnyl);
	geoBuild->createMedium(paraterphnyl);

  // FairGeoMedium* podscintillator = geoMedia->getMedium(esbroot::geometry::superfgd::materials::podscintillator);
	// geoBuild->createMedium(podscintillator);

  // FairGeoMedium* polyethylene = geoMedia->getMedium(esbroot::geometry::superfgd::materials::polyethylene);
	// geoBuild->createMedium(polyethylene);

  // FairGeoMedium* poduleEpoxy = geoMedia->getMedium(esbroot::geometry::superfgd::materials::poduleEpoxy);
	// geoBuild->createMedium(poduleEpoxy);

  // FairGeoMedium* polycarbonate = geoMedia->getMedium(esbroot::geometry::superfgd::materials::polycarbonate);
	// geoBuild->createMedium(polycarbonate);

  // FairGeoMedium* carbonFiber = geoMedia->getMedium(esbroot::geometry::superfgd::materials::carbonFiber);
	// geoBuild->createMedium(carbonFiber);

  FairGeoMedium* fiberCore = geoMedia->getMedium(esbroot::geometry::superfgd::materials::fiberCore);
	geoBuild->createMedium(fiberCore);

  FairGeoMedium* fiberCladding = geoMedia->getMedium(esbroot::geometry::superfgd::materials::fiberCladding);
	geoBuild->createMedium(fiberCladding);

  FairGeoMedium* fairTiO2 = geoMedia->getMedium(esbroot::geometry::superfgd::materials::titaniumDioxide);
  geoBuild->createMedium(fairTiO2);

  FairGeoMedium* fairPolystyrene = geoMedia->getMedium(esbroot::geometry::superfgd::materials::polystyrene);
  geoBuild->createMedium(fairPolystyrene);

  FairGeoMedium* fairAir = geoMedia->getMedium(esbroot::geometry::superfgd::materials::air);
  geoBuild->createMedium(fairAir);

  FairGeoMedium* vacuum = geoMedia->getMedium(esbroot::geometry::superfgd::materials::vacuum);
  geoBuild->createMedium(vacuum);
}

void FgdMCGenFitRecon::PrintFitTrack(genfit::Track& fitTrack)
{
  const genfit::MeasuredStateOnPlane& me = fitTrack.getFittedState();
  LOG(debug)<< "\tFitted Momentum  [" << (me.getMom()).Mag() <<"]" << "(" << (me.getMom()).X() << "," << (me.getMom()).Y() << "," << (me.getMom()).Z() << ")";
  //LOG(debug)<< " X  " << (me.getMom()).X()<< " Y " << (me.getMom()).Y()<< " Z  " << (me.getMom()).Z();

  genfit::FitStatus* fiStatuStatus = fitTrack.getFitStatus();
  fiStatuStatus->Print();

  LOG(debug)<< "fiStatuStatus->isFitted()  " << fiStatuStatus->isFitted();
  LOG(debug)<< "fiStatuStatus->isFitConverged()  " << fiStatuStatus->isFitConverged();
  LOG(debug)<< "fiStatuStatus->isFitConvergedFully()  " << fiStatuStatus->isFitConvergedFully();
  LOG(debug)<< "fiStatuStatus->isFitConvergedPartially()  " << fiStatuStatus->isFitConvergedPartially();
  LOG(debug)<< "fitTrack.getNumPoints() " << fitTrack.getNumPoints();
}

void FgdMCGenFitRecon::WriteOutput( Int_t pdg
                          , const TVector3& fitMom
                          , const TVector3& mcMom
                          , const genfit::Track& fitTrack
                          , genfit::FitStatus*& fiStatuStatus
                          , std::vector<ReconHit>& track)
{
  Float_t&& temp = fitMom.Mag() - mcMom.Mag();
  //LOG(WARNING)<< "Diff " << " Fitted " << fitMom.Mag() << " - Mc " << mcMom.Mag() << " = " << temp;
  if(pdg == genie::kPdgElectron || pdg == genie::kPdgPositron)      { felectronFitErr.emplace_back(temp);}
  if(pdg == genie::kPdgMuon || pdg == genie::kPdgAntiMuon)          { fmuonFitErr.emplace_back(temp); }
  if(pdg == genie::kPdgProton)                                      { fprotonFitErr.emplace_back(temp);}
  if(pdg == genie::kPdgPiP || pdg == genie::kPdgPiM)                { fpionFitErr.emplace_back(temp);}


  std::set<Int_t> xproj;
  std::set<Int_t> yproj;
  std::set<Int_t> zproj;
  for(Int_t i = 0; i < track.size(); ++i)
  {
      xproj.insert(track[i].fmppcLoc.X());
      yproj.insert(track[i].fmppcLoc.Y());
      zproj.insert(track[i].fmppcLoc.Z());
  }
  FitData dataItem;
  dataItem.fitMom = fitMom;
  dataItem.mcMom = mcMom;
  dataItem.trackProj.SetXYZ(xproj.size(), yproj.size(), zproj.size());
  dataItem.totalPath = track.size();
  dataItem.isFitted = fiStatuStatus->isFitted() && fiStatuStatus->isFitConverged();

  if(pdg == genie::kPdgElectron || pdg == genie::kPdgPositron)      { felectronFitData.emplace_back(dataItem);}
  if(pdg == genie::kPdgMuon || pdg == genie::kPdgAntiMuon)          { fmuonFitData.emplace_back(dataItem); }
  if(pdg == genie::kPdgProton)                                      { fprotonFitData.emplace_back(dataItem);}
  if(pdg == genie::kPdgPiP || pdg == genie::kPdgPiM)                { fpionFitData.emplace_back(dataItem);}
}

Bool_t FgdMCGenFitRecon::isParticleNeutral(Int_t pdg)
{
  Bool_t isNeutral =  (pdg ==  genie::kPdgNeutron) ||
                      (pdg ==  genie::kPdgPi0) ||
                      (pdg ==  genie::kPdgGamma) ||
                      genie::pdg::IsNeutralLepton(pdg);

  return isNeutral;
}

Bool_t FgdMCGenFitRecon::isAllowed(Int_t pdg)
{
  Bool_t allowed = (pdg ==  genie::kPdgElectron) ||
                      (pdg ==  genie::kPdgPositron) ||
                      (pdg ==  genie::kPdgMuon) ||
                      (pdg ==  genie::kPdgAntiMuon) ||
                      (pdg ==  genie::kPdgProton) ||
                      (pdg ==  genie::kPdgPiP) ||
                      (pdg ==  genie::kPdgPiM);

  return allowed;
}

Long_t FgdMCGenFitRecon::ArrInd(int x, int y, int z)
{
  return (x*f_bin_Y*f_bin_Z + y*f_bin_Z+z);
}

// -------------------------------------------------------------------------

}// namespace superfgd
}// namespace reconstruction
}// namespace esbroot
