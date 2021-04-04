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
#include "ProlateSpacepointMeasurement.h"
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
#include <cmath>
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
    Float_t electron_threashold_cut = 1.6;
    TTree * electronFitErrTree = new TTree("electronFitErrTree"
                                ,esbroot::geometry::superfgd::DP::FGD_TMVA_DATA_ROOT_FILE.c_str());

    Float_t electronFitErr = 0.;
    electronFitErrTree->Branch("electronFitErr", &electronFitErr);
    for(size_t i = 0; i < felectronFitErr.size(); ++i)
    {
      if(felectronFitErr[i] < electron_threashold_cut)
      {
        electronFitErr = felectronFitErr[i];
        electronFitErrTree->Fill();
      }
      
    }

    TTree * electronFitErrTree_Above_Thres = new TTree("electronFitErrTree_Above_Thres"
                                ,esbroot::geometry::superfgd::DP::FGD_TMVA_DATA_ROOT_FILE.c_str());

    Float_t electronFitErr_Above_Thres = 0.;
    electronFitErrTree_Above_Thres->Branch("electronFitErr_Above_Thres", &electronFitErr_Above_Thres);
    for(size_t i = 0; i < felectronFitErr.size(); ++i)
    {
      if(electron_threashold_cut < felectronFitErr[i])
      {
        electronFitErr_Above_Thres = felectronFitErr[i];
        electronFitErrTree_Above_Thres->Fill();
      }
      
    }
    //========================================================

    //========================================================
    Float_t muon_threashold_cut = 1.6;
    TTree * muonFitErrTree = new TTree("muonFitErrTree"
                                ,esbroot::geometry::superfgd::DP::FGD_TMVA_DATA_ROOT_FILE.c_str());

    Float_t muonFitErr = 0;
    muonFitErrTree->Branch("muonFitErr", &muonFitErr);
    for(size_t i = 0; i < fmuonFitErr.size(); ++i)
    {
      
      if(fmuonFitErr[i] < muon_threashold_cut)
      {
          muonFitErr = fmuonFitErr[i];
          muonFitErrTree->Fill();
      }      
    }

    TTree * muonFitErrTree_Above_Thres = new TTree("muonFitErrTree_Above_Thres"
                                ,esbroot::geometry::superfgd::DP::FGD_TMVA_DATA_ROOT_FILE.c_str());

    Float_t muonFitErrAboveThr = 0;
    muonFitErrTree_Above_Thres->Branch("muonFitErr_Above_Threashold", &muonFitErrAboveThr);
    for(size_t i = 0; i < fmuonFitErr.size(); ++i)
    {
      
      if(muon_threashold_cut < fmuonFitErr[i])
      {
          muonFitErrAboveThr = fmuonFitErr[i];
          muonFitErrTree_Above_Thres->Fill();
      }      
    }
    //========================================================

    //========================================================
    TTree * muonFitErrXYZTree = new TTree("muonFitErrXYZTree"
                                ,esbroot::geometry::superfgd::DP::FGD_TMVA_DATA_ROOT_FILE.c_str());

    Float_t muonFitXErr = 0;
    Float_t muonFitYErr = 0;
    Float_t muonFitZErr = 0;
    muonFitErrXYZTree->Branch("muonFitXErr", &muonFitXErr);
    muonFitErrXYZTree->Branch("muonFitYErr", &muonFitYErr);
    muonFitErrXYZTree->Branch("muonFitZErr", &muonFitZErr);
    for(size_t i = 0; i < fmuonFitData.size(); ++i)
    {
      FitData& dat = fmuonFitData[i];
      muonFitXErr = dat.fitMom.X() - dat.mcMom.X();
      muonFitYErr = dat.fitMom.Y() - dat.mcMom.Y();
      muonFitZErr = dat.fitMom.Z() - dat.mcMom.Z();
      muonFitErrXYZTree->Fill();
      muonFitErrXYZTree->Fill();
      muonFitErrXYZTree->Fill();
    }
    //========================================================

    //========================================================
    Float_t proton_threashold_cut = 1.6;
    TTree * protonFitErrTree = new TTree("protonFitErrTree"
                                ,esbroot::geometry::superfgd::DP::FGD_TMVA_DATA_ROOT_FILE.c_str());

    Float_t protonFitErr = 0;;
    protonFitErrTree->Branch("protonFitErr", &protonFitErr);
    for(size_t i = 0; i < fprotonFitErr.size(); ++i)
    {
      if(fprotonFitErr[i] < proton_threashold_cut)
      {
        protonFitErr = fprotonFitErr[i];
        protonFitErrTree->Fill();
      }
    }

    TTree * protonFitErrTree_Above_Threashold = new TTree("protonFitErrTree_Above_Threashold"
                                ,esbroot::geometry::superfgd::DP::FGD_TMVA_DATA_ROOT_FILE.c_str());

    Float_t protonFitErr_Above_Threashold = 0;;
    protonFitErrTree_Above_Threashold->Branch("protonFitErr_Above_Threashold", &protonFitErr_Above_Threashold);
    for(size_t i = 0; i < fprotonFitErr.size(); ++i)
    {
      if(proton_threashold_cut < fprotonFitErr[i])
      {
        protonFitErr_Above_Threashold = fprotonFitErr[i];
        protonFitErrTree_Above_Threashold->Fill();
      }
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
    outFile->WriteTObject(electronFitErrTree_Above_Thres);

    outFile->WriteTObject(muonFitErrTree);
    outFile->WriteTObject(muonFitErrTree_Above_Thres);
    outFile->WriteTObject(muonFitErrXYZTree);  

    outFile->WriteTObject(protonFitErrTree); 
    outFile->WriteTObject(protonFitErrTree_Above_Threashold); 
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

    writeErrFile(std::string("_muon"),        fmuonFitData);
    writeErrFile(std::string("_proton"),      fprotonFitData);
    writeErrFile(std::string("_electron"),    felectronFitData);
    writeErrFile(std::string("_pion"),        fpionFitData);

    writeTrErrFile(std::string("_muon"),        fmuonFitData);
    writeTrErrFile(std::string("_proton"),      fprotonFitData);
    writeTrErrFile(std::string("_electron"),    felectronFitData);
    writeTrErrFile(std::string("_pion"),        fpionFitData);
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

  // std::map<int, std::vector<float>> projErr;
  // for(size_t i = 0; i < dataVec.size(); ++i)
  // {
  //   FitData& data = dataVec[i];
  //   if(!data.isFitted) continue;

  //   int proj(0);
  //   float err(0.);
  //   switch(axis)
  //   {
  //     case Axis::X: proj = data.trackProj.X();
  //                   err = (data.mcMom.X() - data.fitMom.X());
  //                   break;
  //     case Axis::Y: proj = data.trackProj.Y();
  //                   err = (data.mcMom.Y() - data.fitMom.Y());
  //                   break;
  //     case Axis::Z: proj = data.trackProj.Z(); 
  //                   err = (data.mcMom.Z() - data.fitMom.Z());;
  //                   break;
  //     default: 
  //                   break;
  //   }

  //   bool hasProj = (projErr.find(proj) != projErr.end());
  //   if(hasProj)
  //   {
  //     projErr[proj].insert(projErr[proj].end(), err);
  //   }
  //   else
  //   {
  //     std::vector<float> errs;
  //     errs.insert(errs.end(), err);
  //     projErr[proj] = std::move(errs);
  //   }
  // }

  std::ofstream outfile(path, std::ios::trunc);
  if(outfile.is_open())
  {
    for(size_t i = 0; i < dataVec.size(); ++i)
    {
      FitData& data = dataVec[i];
      if(!data.isFitted) continue;

      switch(axis)
      {
        case Axis::X: outfile << data.mcMom.X() << " " << (data.mcMom.X() - data.fitMom.X()) << endl;
                      break;
        case Axis::Y: outfile << data.mcMom.Y() << " " << (data.mcMom.Y() - data.fitMom.Y()) << endl;
                      break;
        case Axis::Z: outfile << data.mcMom.Z() << " " << (data.mcMom.Z() - data.fitMom.Z()) << endl;
                      break;
        default: 
                      break;
      }
    }
    // auto it = projErr.begin();
    // while(it!= projErr.end())
    // {
    //   int keyProj = it->first;
    //   float mean(0.);
    //   float stdDev = calcStdDev(mean, it->second);
    //   ++it;

    //   outfile << keyProj << " " << mean << " " << stdDev << endl;
    // }
  }
  outfile.close();
}

void   FgdMCGenFitRecon::writeErrFile(const std::string& fileEnding, std::vector<FitData>& dataVec)
{
  std::stringstream ss;
  ss << ferrOutPath << fileEnding << ".txt";
  std::string path = ss.str();
  std::ofstream outfile(path, std::ios::trunc);

  if(outfile.is_open())
  {
    for(size_t i = 0; i < dataVec.size(); ++i)
    {
      FitData& data = dataVec[i];
      if(!data.isFitted) continue;

      outfile << data.mcMom.Mag() << " " << (data.mcMom.Mag() - data.fitMom.Mag()) << endl;
    }
  }
  outfile.close(); 
}

void   FgdMCGenFitRecon::writeTrErrFile(const std::string& fileEnding, std::vector<FitData>& dataVec)
{
  std::stringstream ss;
  ss << ferrOutPath << fileEnding << "tr.txt";
  std::string path = ss.str();
  std::ofstream outfile(path, std::ios::trunc);

  if(outfile.is_open())
  {
    for(size_t i = 0; i < dataVec.size(); ++i)
    {
      FitData& data = dataVec[i];
      if(!data.isFitted) continue;

      outfile << data.totalPath << " " << (data.mcMom.Mag() - data.fitMom.Mag()) << endl;
    }
  }
  outfile.close(); 
}

float  FgdMCGenFitRecon::calcStdDev(float& mean, std::vector<float>& vals)
{
  float sum(0.);
  for(size_t i = 0; i < vals.size(); ++i)
  {
    sum+= vals[i];
  }
  mean = sum/vals.size();

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



      // double sum = 0.;
      // for(Int_t bh = 0; bh < hitsOnTrack.size(); ++bh)
      // {
      //   std::sort(hitsOnTrack.begin(), hitsOnTrack.end(), [](ReconHit& bh1, ReconHit& bh2){return bh1.ftime<bh2.ftime;});
      //   ReconHit& hit = hitsOnTrack[bh];
      //   sum += hit.fEdep;
      //   LOG(warning) << " \tMom [ "  << hit.fpdg << " ] \t" << hit.fmom.Mag()  << "  \tdEdx " << hit.fEdep
      //               << "  \tL " << hit.ftrackLength << "  \tE/L " << hit.fEdep/hit.ftrackLength
      //               << "(" << hit.fmppcLoc.X() << "," << hit.fmppcLoc.Y() << "," << hit.fmppcLoc.Z() << ")";
      // }

      //LOG(warning) << " \tAvg dEdx " << 12.71 * sum/hitsOnTrack[hitsOnTrack.size()-1].ftrackLengthOrigin;

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

      // genfit::MaterialEffects::getInstance()->setEnergyLossBetheBloch(false);
      // genfit::MaterialEffects::getInstance()->setMscModel("Highland");
      // genfit::MaterialEffects::getInstance()->drawdEdx(pdg);
      // genfit::MaterialEffects::getInstance()->setNoEffects();
      // genfit::MaterialEffects::getInstance()->setNoiseBetheBloch(false);
      // genfit::MaterialEffects::getInstance()->setNoiseCoulomb(false);
      // genfit::MaterialEffects::getInstance()->setEnergyLossBrems(false);
      // genfit::MaterialEffects::getInstance()->setNoiseBrems(false);
      // genfit::MaterialEffects::getInstance()->ignoreBoundariesBetweenEqualMaterials(false);


      // if(isParticleNeutral(pdg))
      // {
      //   LOG(debug) << "Track " << i << " is of neutral particle ["<< pdg << "] continue with next track.";
      //   continue;
      // }

      // approximate covariance
      static const double resolution = 1;// Default in example is 0.1;   1/sqrt(12)
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
  
      genfit::Track* toFitTrack = new genfit::Track(rep, seedState, seedCov);
      //genfit::Track* toFitTrack = new genfit::Track(rep, posM, momM);

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

        // TMatrixDSym phhitCov(3);
        // // The more photons, the more likely the particle went through the center
        // // dimention of resolution is in [cm]
        // double cube_diag_length = (sqrt(3)/2);
        // double normCoeff = cube_diag_length*(1 - sumPhotons/maxPhotonsPerCube);     //std::sqrt(sumPhotons);
        // //normCoeff = std::sqrt(normCoeff);
        // normCoeff = (normCoeff == 0) ? 1 : normCoeff;
        // phhitCov(0,0) = resolution * resolution  * normCoeff;
        // phhitCov(1,1) = resolution * resolution  * normCoeff;
        // phhitCov(2,2) = resolution * resolution  * normCoeff;


        genfit::AbsMeasurement* measurement = new genfit::SpacepointMeasurement(hitPos, hitCov, detId, bh, nullptr);
        //genfit::AbsMeasurement* measurement = new genfit::ProlateSpacepointMeasurement(hitPos, hitCov, detId, bh, nullptr);
        
        //genfit::AbsMeasurement* measurement = new genfit::SpacepointMeasurement(hitPos, phhitCov, detId, bh, nullptr);
        //measurements.emplace_back(measurement);

        //genfit::AbsMeasurement* measurement = new genfit::SpacepointMeasurement(hitPos, phhitCov, detId, bh, nullptr);
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


        // Double_t momentum=0;
        // bool houghFlag = GetHoughMomentum(hitsOnTrack,momentum);
        // //if()
        // {
        //   float mcMom = std::sqrt(  hitsOnTrack[0].fmom.X()*hitsOnTrack[0].fmom.X()  +  hitsOnTrack[0].fmom.Z()*hitsOnTrack[0].fmom.Z()   );
        //   LOG(debug) << "Found Hough " << houghFlag <<" Hough Momentum " << momentum << " MC mom " << mcMom;
        // }

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
  // float valLimit = 400;
  // if( (pdg == genie::kPdgMuon || pdg == genie::kPdgAntiMuon) && (fitMom.X() > valLimit
  //   || fitMom.Y() > valLimit
  //   || fitMom.Z() > valLimit
  //   || mcMom.X() > valLimit
  //   || mcMom.Y() > valLimit
  //   || mcMom.Z() > valLimit)
  //   )
  //   {
  //     LOG(ERROR) << "Too big momentum?";
  //     exit(1);
  //   }

  bool isFitted = fiStatuStatus->isFitted() && fiStatuStatus->isFitConverged();
  // bool momentumToLowX(false);
  // bool momentumToLowY(false);
  // bool momentumToLowZ(false);
  // static const TVector3 magField = fgdConstructor.GetMagneticField();
  // static const float lowlimit = 0.08;
  // if(magField.X()!=0){ momentumToLowX = (std::abs(mcMom.X())<lowlimit); }
  // if(magField.Y()!=0){ momentumToLowY = (std::abs(mcMom.Y())<lowlimit); }
  // if(magField.Z()!=0){ momentumToLowZ = (std::abs(mcMom.Z())<lowlimit); }
  // bool momentumLow = momentumToLowX || momentumToLowY || momentumToLowZ;
  

  if(!isFitted)
  {
    return;
  }

  Float_t&& temp = fitMom.Mag() - mcMom.Mag();
  // if((pdg == genie::kPdgMuon || pdg == genie::kPdgAntiMuon)    && (temp>4.0 || temp<-4.0))
  // {
  //   LOG(warning)<< "momentumToLowX " << momentumToLowX <<  " momentumToLowX " << momentumToLowY << " momentumToLowZ " << momentumToLowZ;
  //   LOG(warning)<< "std::abs(mcMom.Y()) " << std::abs(mcMom.Y());
  //   LOG(fatal)<< "Diff " << " Fitted " << fitMom.Mag() << " - Mc " << mcMom.Mag() << " = " << temp;
  //   exit(1);
  // }
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


Long_t FgdMCGenFitRecon::HitId(ReconHit& hit)
{
  hit.fHitId = (hit.fmppcLoc.X()*f_bin_Y*f_bin_Z + hit.fmppcLoc.Y()*f_bin_Z + hit.fmppcLoc.Z());
  return hit.fHitId;
}

void FgdMCGenFitRecon::FillHitIds(std::vector<ReconHit>& track )
{
    for(size_t i = 0; i < track.size(); ++i)
    {
      HitId(track[i]);
    }
} 

// Calculate the radius from 3 points using the "Menger curvature" theorem
Double_t FgdMCGenFitRecon::Radius(const TVector3& p1, const TVector3& p2, const TVector3& p3)
{
  //
  //          p1
  //          /\\
  //     y   /  \\  x
  //        /    \\
  //       /______\\
  //      p2        p3
  //          z
  //
  //                    2 sin (Angle(xyz))
  //      c(x,y,z) = ------------------------
  //                        |x - z|
  //
  //    where Angle(xyz) is the angle from point p2 to point p1 and p3
  //    and |x - z| is the magnitude of the line connecting points p2 and p3

  TVector3 x = p3 - p1;
  TVector3 y = p2 - p1;
  TVector3 z = p2 - p3;

  Double_t angle = y.Angle(x);
  Double_t x_z_Mag = z.Mag();

  Double_t c = (2*std::sin(angle))/x_z_Mag;
  Double_t R = 1./c;

  LOG(debug3) << "Radius is -> " << R << " [cm]";

  return R;
}

bool FgdMCGenFitRecon::GetHoughMomentum(std::vector<ReconHit>& track, Double_t& momentum)
{
    // For calculation charge is taken as 1 unit of 'e'
    static const Double_t charge = 0.3; // Charge is 1e, but to convert from meter * T to GeV/c a factor of 0.3 has to be added
                                        // thus 1e * 0.3 => 0.3
    static const Double_t inf = std::numeric_limits<Double_t>::infinity();

    bool rc(false);
    TVector3 p1 = track[0].fmppcLoc;
    TVector3 p2 = track[track.size()/2].fmppcLoc;
    TVector3 p3 = track[track.size()-1].fmppcLoc;

    TVector3 magField = fgdConstructor.GetMagneticField();

    size_t countAxis{0};

    // We don`t want to calculate the momentum in the direction of the magnetic field
    if(magField.X()!=0)
    {
      p1.SetX(0);
      p2.SetX(0);
      p3.SetX(0);
      ++countAxis;
    }

    // We don`t want to calculate the momentum in the direction of the magnetic field
    if(magField.Y()!=0)
    {
      p1.SetY(0);
      p2.SetY(0);
      p3.SetY(0);
      ++countAxis;
    }

    // We don`t want to calculate the momentum in the direction of the magnetic field
    if(magField.Z()!=0)
    {
      p1.SetZ(0);
      p2.SetZ(0);
      p3.SetZ(0);
      ++countAxis;
    }

    // We can calculate the momentum if we have a magnetic field only in ONE axis
    // the radius is extracted from a 2D plane, not a 3D. 
    if(countAxis!=1)
    {
      return false;
    }

    Double_t radius = Radius(p1, p2, p3); // get initial radius
    Double_t houghRadius = 0;

    Double_t temp = 0;
    Int_t cubesInCirlcle = 0;
    if(!std::isnan(radius) 
        && radius!=inf 
        && HoughRadius(track, temp, houghRadius, cubesInCirlcle))
    // Double_t temp = 0;
    // if(HoughRadius(track, temp, houghRadius))
    {
      rc = true;

      // The magnetic field is expected to be only in one axis
      Double_t magField_T =  magField.X();
      if(magField.Y()!=0) {magField_T = magField.Y();}
      if(magField.Z()!=0) {magField_T = magField.Z();}
      magField_T = magField_T / 10.; // convert from kGauss to Tesla units

      Double_t totalEdep = 0;
      for(int i = 0; i < track.size(); i+=2)
      {
        totalEdep += track[i].fEdep;
      }


      Double_t R = houghRadius/100.; // convert in meters
      
      momentum = charge * R * magField_T;
      LOG(debug) << "Menger [r] = " << radius << ", Hough [r] = " << houghRadius;
      LOG(debug) << "Menger [p] = " << (charge * (radius/100) * magField_T) << ", Hough [p] = " 
                                    << (charge * R * magField_T) 
                                    << " => total Edep = " << totalEdep;

      // float partCubesCovered = (float)cubesInCirlcle/track.size();
      // LOG(debug) << "Menger [p] = " << (charge * (radius/100) * magField_T) << ", Hough [p] = " 
      //                               << (partCubesCovered * charge * R * magField_T) 
      //                               << " => total Edep = " << totalEdep;
    }

    return rc;
}

bool FgdMCGenFitRecon::HoughRadius(std::vector<ReconHit>& track, Double_t& initialRadius, Double_t& radius, Int_t& cubesInCirlcle)
{
  // Hough circle transform
  // For each A[a,b,r] = 0; // fill with zeroes initially, instantiate 3D matrix
  // For each cell(x,y)
  //     For each theta t = 0 to 360  // the possible  theta 0 to 360
  //        b = y – r * sin(t * PI / 180);  //polar coordinate for center (convert to radians)
  //        a = x – r * cos(t * PI / 180); //polar coordinate for center (convert to radians)
  //        A[a,b,r] +=1; //voting
  //     end
  // end

  // Double_t rot90deg = TMath::Pi()/2;

  bool rc(false);

  FillHitIds(track);
  TVector3 magField = fgdConstructor.GetMagneticField();
  auto degToRad = (TMath::Pi() / 180);
  std::map<int, std::map<long long, HoughPoint>> houghPoints;

  Quadrant quad = Quadrant::One;
  SemiCircle semi = SemiCircle::SemiOne;
  // Find quadrant and semicircle
  {
    ReconHit& first = track[0];
    ReconHit& last = track[track.size() -1];
    int x = 0;
    int y = 0;
    if(magField.X()!=0)
    {
        x = last.fmppcLoc.Y() - first.fmppcLoc.Y();
        y = last.fmppcLoc.Z() - first.fmppcLoc.Z();
    }
    else if(magField.Y()!=0)
    {
        x = last.fmppcLoc.X() - first.fmppcLoc.X();
        y = last.fmppcLoc.Z() - first.fmppcLoc.Z();
    }
    else
    {
        x = last.fmppcLoc.X() - first.fmppcLoc.X();
        y = last.fmppcLoc.Y() - first.fmppcLoc.Y();
    }

    if(y >= 0) {semi = SemiCircle::SemiOne;}
    if(y < 0) {semi = SemiCircle::SemiTwo;}

    if(y >= 0) 
    { 
        if(x >= 0) {quad = Quadrant::One;}
        else{quad = Quadrant::Two;}
    }
    if(y < 0) 
    { 
        if(x >= 0) {quad = Quadrant::Four;}
        else{quad = Quadrant::Three;}
    }
  }
  
  int cubeLength = (f_step_X < f_step_Y) ? f_step_Y : f_step_X;
  cubeLength = (cubeLength < f_step_Z) ? f_step_Z : cubeLength;
  
  Double_t range = 500;
  int min_r = cubeLength * 1;

  int r_start = initialRadius - range;
  r_start = (r_start < min_r) ? min_r : r_start;
  int r_end = initialRadius + range;
  int r_delta = 1;

  int t_delta = 1;

  int t_start = 0;
  int t_end = 360;

  switch(semi)
  { 
    case SemiCircle::SemiOne: t_start = 1; t_end = 179; 
                          break;
    case SemiCircle::SemiTwo: t_start = 181; t_end = 359;
                          break;
    default:
          break;
  }


  for(int r = r_start;  r <= r_end; r+= r_delta)
  {
    std::map<long long, HoughPoint> currHoughPoints;

    for(int i = 0; i < (track.size() * 0.8); i+=1)
    {
      int x = 0;
      int y = 0;
      // Assume magnetic field is only in one axis
      ReconHit& hit = track[i];
      if(magField.X()!=0)
      {
          x = hit.fmppcLoc.Y();
          y = hit.fmppcLoc.Z();
      }
      else if(magField.Y()!=0)
      {
          x = hit.fmppcLoc.X();
          y = hit.fmppcLoc.Z();
      }
      else
      {
          x = hit.fmppcLoc.X();
          y = hit.fmppcLoc.Y();
      }
      for(int t = t_start;  t <= t_end; t+=t_delta)
      {
          HoughPoint point;
          point.a = x - r * cos(t * degToRad);  //polar coordinate for center (convert to radians)
          point.b = y - r * sin(t * degToRad);  //polar coordinate for center (convert to radians)
          point.r = r;
          
          long long id = (long long)point.a*r_end*r_end + (long long)point.b*r_end + point.r;

          std::map<long long, HoughPoint>::iterator it = currHoughPoints.find(id);
          if(it == currHoughPoints.end())
          {
            currHoughPoints[id] = point;
            point.count.insert(hit.fHitId);
          }else
          {
            HoughPoint& p = currHoughPoints[id];
            p.count.insert(hit.fHitId);
          }
      }
    }
    houghPoints[r] = currHoughPoints;
  }

  int currentMax = -1;
  HoughPoint maxpoint;
  std::map<int, std::map<long long, HoughPoint>>::iterator it = houghPoints.begin();
  while( it!= houghPoints.end())
  {
      
      std::map<long long, HoughPoint> tempMap = it->second;
      // LOG(debug) << "Radius "<<  it->first << "; Entries " << tempMap.size();
       
      std::map<long long, HoughPoint>::iterator tempIt = tempMap.begin();
      while( tempIt!= tempMap.end())
      {
        HoughPoint& point = tempIt->second;
        // LOG(debug) << "point.count "<<  point.count.size();
        const int& sizeCount = point.count.size();
        if( currentMax < sizeCount)
        {
            currentMax = sizeCount;
            maxpoint = point;
        }
        ++tempIt;
      }

      ++it;
  }

  it = houghPoints.begin();
  while( it!= houghPoints.end())
  {
      
      std::map<long long, HoughPoint> tempMap = it->second;
      // LOG(debug) << "Radius "<<  it->first << "; Entries " << tempMap.size();
       
      std::map<long long, HoughPoint>::iterator tempIt = tempMap.begin();
      while( tempIt!= tempMap.end())
      {
        HoughPoint& point = tempIt->second;
        // LOG(debug) << "point.count "<<  point.count.size();
        const int& sizeCount = point.count.size();
        if( currentMax == sizeCount && maxpoint.r > point.r)
        {
            maxpoint = point;
        }
        ++tempIt;
      }

      ++it;
  }

  if(currentMax != -1)
  {
    rc = true;
    radius = maxpoint.r;
    cubesInCirlcle = maxpoint.count.size();
    LOG(debug) << "Found A["<< maxpoint.a << " ," << maxpoint.b << " ," << maxpoint.r << "] = " << maxpoint.count.size();
  }

  return rc;
}
// -------------------------------------------------------------------------

}// namespace superfgd
}// namespace reconstruction
}// namespace esbroot
