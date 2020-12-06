#include "EsbReconstruction/EsbSuperFGD/FgdGraphGenFitRecon.h"
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
#include <TFile.h>

//Genie headers
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
FgdGraphGenFitRecon::FgdGraphGenFitRecon() :
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
  , feventNum(0)
  , fuseSmoothPos(false)
{ 
}
// -------------------------------------------------------------------------

// -----   Constructor   -------------------------------------------
FgdGraphGenFitRecon::FgdGraphGenFitRecon(const char* name
                          , const char* geoConfigFile
                          , const char* graphConfig
                          , const char* mediaFile
                          , const char* eventData
                          , const char* outputRootFile
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
  , fgeoConfFile(geoConfigFile)
  , fgraphConfig(graphConfig)
  , fmediaFile(mediaFile)
  , feventData(eventData)
  , foutputRootFile(outputRootFile)
  , fTracksArray(nullptr)
  , fdisplay(nullptr)
  , isGenFitVisualization(visualize)
  , fGenFitVisOption(visOption)
  , fminGenFitInterations(2)
  , fmaxGenFitIterations(4)
  , fminHits(25)
  , feventNum(0)
  , fuseSmoothPos(false)
  , freconTemplate(geoConfigFile, graphConfig)
{ 
  fParams.LoadPartParams(geoConfigFile);
}
// -------------------------------------------------------------------------



// -----   Destructor   ----------------------------------------------------
FgdGraphGenFitRecon::~FgdGraphGenFitRecon() 
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
InitStatus FgdGraphGenFitRecon::Init() 
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

  // init geometry and mag. field
  TVector3 magField = fgdConstructor.GetMagneticField(); // values are in kGauss
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(magField.X(),magField.Y(), magField.Z())); 
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
  genfit::MaterialEffects::getInstance()->setEnergyLossBetheBloch(true);
  genfit::MaterialEffects::getInstance()->setEnergyLossBrems(true);
  genfit::MaterialEffects::getInstance()->setDebugLvl(fDebuglvl_genfit);

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

  std::ifstream eventFileStream;
  try
  {        
      eventFileStream.open(feventData.c_str(), std::ios::in);

      if(eventFileStream.is_open())
      {
          std::string line;
          while(std::getline(eventFileStream,line))
          {
            feventRecords.emplace_back(FgdTMVAEventRecord(line));
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
  

  return kSUCCESS;
}

void FgdGraphGenFitRecon::OutputFileInit(FairRootManager* manager)
{
  // Create and register output array
  fTracksArray = new TClonesArray(genfit::Track::Class(), 1000);
  manager->Register(geometry::superfgd::DP::FGD_FIT_TRACK.c_str()
                    , geometry::superfgd::DP::FGD_BRANCH_FIT.c_str()
                    , fTracksArray, kTRUE);
}


// -------------------------------------------------------------------------



// -----   Public methods   --------------------------------------------
void FgdGraphGenFitRecon::FinishEvent()
{
  
}

void FgdGraphGenFitRecon::FinishTask()
{
  if(isGenFitVisualization)
  {
    fdisplay->open();
  }

  static const TVector3 z_axis(0,0,1);

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
    //========================================================


    const Int_t evInd = feventRecords.size();
    
    Int_t limit = evInd < feventNum ? evInd : feventNum;

    FgdTMVAEventRecord* dataEvent = nullptr;
    for(size_t ind = 0 ; ind < limit; ind++)
    {
        LOG(debug2) << "Writing data for event " << ind;
        dataEvent = &feventRecords[ind];

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
        fit_muon_mom =  dataEvent->GetMuonFitMom().Mag();
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


    outFile->WriteTObject(trainTree);  
    outFile->WriteTObject(fittedMomTree);   
    outFile->WriteTObject(calorimetricMomTree);                   
  }
  catch(...)
  {

  }
  
  outFile->Close();
  
  delete outFile;
}

void FgdGraphGenFitRecon::Exec(Option_t* opt) 
{  
  try
  {

    FgdTMVAEventRecord& tvmaEventRecord = feventRecords[feventNum];
    tvmaEventRecord.ReadEventData();

    std::vector<ReconHit> allhits;
    std::vector<std::vector<ReconHit>> foundTracks;

    bool rc = GetHits(allhits);
    if(rc)
    { 
      LOG(debug) <<" Hits to retrieve tracks from " << allhits.size();
      rc = FindUsingGraph(allhits, foundTracks);
    }

    if(rc)
    {
      LOG(debug) <<" Tracks found " << foundTracks.size();
      std::vector<ReconHit>& longestTrack = foundTracks[0];
      size_t ltId(0);
      for(size_t i = 0; i <  foundTracks.size() ; ++i)
      {
        std::vector<ReconHit>& track = foundTracks[i];
        std::sort(track.begin(), track.end(), [](ReconHit& bh1, ReconHit& bh2){return bh1.fmppcLoc.Z()<bh2.fmppcLoc.Z();});
        LOG(debug) << "size " << track.size();
        LOG(debug) << "start \tX " << track[0].fmppcLoc.X() << " Y " << track[0].fmppcLoc.Y()<< " Z " << track[0].fmppcLoc.Z();
        LOG(debug) << "end \tX " << track[track.size()-1].fmppcLoc.X() << " Y " << track[track.size()-1].fmppcLoc.Y()<< " Z " << track[track.size()-1].fmppcLoc.Z();
        LOG(debug) << "=====";
      }
      FitTracks(foundTracks);
    }
    else
    {
      LOG(error) << " Could not find clean hits or tracks! ";
    }

    if(rc)
    {
      LOG(debug) <<" Extracting tmva data ";
      ExtractTMVAdata(allhits);
    }
    else
    {
      LOG(error) << " Could not extract tmva data ";
    }
  }
  catch(genfit::Exception& e)
  {
      LOG(fatal) << "Exception, when tryng to fit track";
      LOG(fatal) << e.what();
  }
  ++feventNum;
}
// -------------------------------------------------------------------------


// -----   Protected methods   --------------------------------------------

Bool_t FgdGraphGenFitRecon::GetHits(std::vector<ReconHit>& allHits)
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

    LOG(debug2) << "TrackId " << hit->GetTrackId();
    LOG(debug2) << "GetPgd " << hit->GetPgd();
    LOG(debug2) << "GetTrackLengthOrigin " << hit->GetTrackLengthOrigin();
    LOG(debug2) << "GetTrackLenght " << hit->GetTrackLenght();
    
    Int_t ind = ArrInd(x,y,z);
    if(visited[ind])
    {
      // If already exists, add the photons
      // ReconHit toFind;
      // toFind.fmppcLoc = mppcLoc;
      // std::vector<ReconHit>::iterator recHit = find(allHits.begin(), allHits.end(), toFind);
      // ReconHit& foundHit = *recHit;
      // foundHit.fphotons = foundHit.fphotons + photoE;
      continue;
    }
    visited[ind] = true;

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

Bool_t FgdGraphGenFitRecon::FindUsingGraph(std::vector<ReconHit>& hits
                  , std::vector<std::vector<ReconHit>>& foundTracks)
{
  if(hits.empty())
  {
    return false;
  }

  freconTemplate.BuildGraph(hits);
  freconTemplate.FindLeaves(hits);

  if(fuseSmoothPos) freconTemplate.SmoothGraph(hits);
  //BuildGraph(hits);

  // Print out the build graph
  for(Int_t i=0; i<hits.size(); ++i)
  {
    LOG(debug2) << "i " << i;
    ReconHit* center = &hits[i];
    for(Int_t j=0; j<hits[i].fAllHits.size(); ++j)
    {
      ReconHit* localHit = hits[i].fAllHits[j];
      TVector3 diff = center->fmppcLoc - localHit->fmppcLoc;
      LOG(debug2) << " Local Id "<< hits[i].fAllHits[j]->fLocalId << " \t X " << diff.X() << " \t Y " << diff.Y() << " \t Z " << diff.Z();
    }
    LOG(debug2) << "X " << hits[i].fmppcLoc.X() << " Y " << hits[i].fmppcLoc.Y()<< " Z " << hits[i].fmppcLoc.Z();
    LOG(debug2) << "Photons "<< " X " << hits[i].fphotons.X() << " Y " << hits[i].fphotons.Y()<< " Z " << hits[i].fphotons.Z();
    LOG(debug2) << "=====";
  }

  std::vector<std::vector<ReconHit*>> tracks;
  freconTemplate.FindTracks(hits,tracks);

  LOG(debug) <<"Leaves found " << tracks.size(); // Initially leaves are equal to the number of tracks

  freconTemplate.CalculateGrad(tracks, fuseSmoothPos);
  std::vector<std::vector<ReconHit*>> splitTracks;
  freconTemplate.SplitTrack(tracks, splitTracks);
  //std::vector<std::vector<ReconHit*>>& splitTracks = tracks;

  Int_t totalHitsInTracks(0);
  for(size_t i = 0; i<splitTracks.size(); ++i)
  {
    std::vector<ReconHit*>& track = splitTracks[i];
    std::vector<ReconHit> currentTrack;
    LOG(debug2) << "=============================== ";
    LOG(debug2) << "Track " << i;
    
    totalHitsInTracks+=track.size();

    for(size_t j = 0; j<track.size(); ++j)
    {
      ReconHit* trackHit = track[j];

      //if(j==0) // first hit
      // {
      //   LOG(debug2) << "Pdg code particle " << trackHit->fpdg;
      //   LOG(debug2) << "Particle momentum " << trackHit->fmom.Mag();
      //   LOG(debug2) << " \tmomX " << trackHit->fmom.X() << " \tmomY " << trackHit->fmom.Y() << " \tmomZ " << trackHit->fmom.Z();
      // }

      currentTrack.emplace_back(*trackHit);   

      LOG(debug2) << " \tX " << trackHit->fmppcLoc.X() << " \tY " << trackHit->fmppcLoc.Y() << " \tZ " << trackHit->fmppcLoc.Z();    

      // if(j== (track.size() -1 )) // last hit
      // {
      //   LOG(debug2) << "Pdg code particle " << trackHit->fpdg;
      //   LOG(debug2) << "Particle momentum " << trackHit->fmom.Mag();
      //   LOG(debug2) << " \tmomX " << trackHit->fmom.X() << " \tmomY " << trackHit->fmom.Y() << " \tmomZ " << trackHit->fmom.Z();
      // }
    }
    LOG(debug2) << "=============================== ";

    foundTracks.push_back(currentTrack);
  }

  LOG(debug) << "Total hits in tracks " << totalHitsInTracks;

  return !foundTracks.empty();
}

// Bool_t FgdGraphGenFitRecon::CalculateInitialMomentum(const std::vector<ReconHit>& track,const TVector3& magField, TVector3& momentum , TVector3& momentumLoss)
// {
//   if(track.empty())
//   {
//       momentum.SetXYZ(0.,0.,0);
//       momentumLoss.SetXYZ(0.,0.,0);
//       return false;
//   }

//   Int_t segment = fParams.ParamAsInt(esbroot::geometry::superfgd::DP::FGD_TRACK_MOMENTUM_SEGMENT);
//   Int_t avgInitialPoints = fParams.ParamAsInt(esbroot::geometry::superfgd::DP::FGD_INITIAL_TRACK_POINTS_MOMENTUM);
//   static const Double_t inf = std::numeric_limits<Double_t>::infinity();

//   const Int_t defaultSegment = 3;
//   segment = (segment<=0) ? defaultSegment : segment;
//   avgInitialPoints = (avgInitialPoints<=0) ? defaultSegment : avgInitialPoints;

//   // If the track  lenght is less than the set segment, calculate the momentum from the arc of the whole track
//   // using for the 3 points the begin, end and all the middle points between them
//   Int_t lengthSize = (track.size() < segment ) ? (track.size() -1) : (segment -1) ;

//   std::vector<TVector3> trackMomentums;
//   std::vector<TVector3> tarckMomentumLosses;

//   TVector3 trackVector = track[track.size()-1].fHitPos - track[0].fHitPos;

//   for(size_t i = lengthSize; i < track.size(); i+=3)
//   {
//     size_t p1_pos = i - lengthSize;
//     TVector3 p1 = track[p1_pos].fHitPos;
//     TVector3 p3 = track[i].fHitPos;

//     size_t p2_pos = p1_pos + 1;
//     //size_t p2_pos = p1_pos + lengthSize/2;
//     std::vector<TVector3> segmentMomentums;

//     while(p2_pos<i)
//     {
//       TVector3 segmentMom(0,0,0);
//       TVector3 p2 = track[p2_pos].fHitPos;
//       if(CalculateMomentum(trackVector, p1, p2, p3, magField, segmentMom)  
//           &&  segmentMom.Mag()!=inf
//           && !std::isnan(momentum.Mag()))
//       {
//         segmentMomentums.push_back(segmentMom);
//       }
//       ++p2_pos;
//     }

//     Double_t xMom(0);
//     Double_t yMom(0);
//     Double_t zMom(0);
//     for(size_t j = 0; j< segmentMomentums.size() ; ++j)
//     {
//       xMom+=segmentMomentums[j].X();
//       yMom+=segmentMomentums[j].Y();
//       zMom+=segmentMomentums[j].Z();
//     }

//     xMom = xMom/segmentMomentums.size();
//     yMom = yMom/segmentMomentums.size();
//     zMom = zMom/segmentMomentums.size();

//     trackMomentums.emplace_back(xMom, yMom, zMom);
//   }

//   for(size_t i = 1; i< trackMomentums.size() ; ++i)
//   {
//     TVector3 momLoss = trackMomentums[i] - trackMomentums[i-1];
//     tarckMomentumLosses.push_back(momLoss);
//   }

//   Double_t xTrackMom(0);
//   Double_t yTrackMom(0);
//   Double_t zTrackMom(0);
//   Int_t limit = avgInitialPoints < trackMomentums.size()? avgInitialPoints : trackMomentums.size();

//   for(size_t i = 0; i < limit ; ++i)
//   {
//     xTrackMom+=trackMomentums[i].X();
//     yTrackMom+=trackMomentums[i].Y();
//     zTrackMom+=trackMomentums[i].Z();;
//   }

//   xTrackMom = xTrackMom/ limit;
//   yTrackMom = yTrackMom/ limit;
//   zTrackMom = zTrackMom/ limit;

//   momentum.SetXYZ(xTrackMom, yTrackMom, zTrackMom);


//   Double_t xTrackMomLoss(0);
//   Double_t yTrackMomLoss(0);
//   Double_t zTrackMomLoss(0);
//   Int_t limitLoss = avgInitialPoints < tarckMomentumLosses.size()? avgInitialPoints : tarckMomentumLosses.size();
//   for(size_t i = 0; i < limitLoss ; ++i)
//   {
//     xTrackMomLoss+=tarckMomentumLosses[i].X();
//     yTrackMomLoss+=tarckMomentumLosses[i].Y();
//     zTrackMomLoss+=tarckMomentumLosses[i].Z();;
//   }

//   xTrackMomLoss = xTrackMomLoss/ limitLoss;
//   yTrackMomLoss = yTrackMomLoss/ limitLoss;
//   zTrackMomLoss = zTrackMomLoss/ limitLoss;

//   momentumLoss.SetXYZ(xTrackMomLoss, yTrackMomLoss, zTrackMomLoss);


//   LOG(debug3) << " =========================== Track Momentum =========================== ";
//   for(size_t i = 0; i < trackMomentums.size() ; ++i)
//   {
//     LOG(debug3) << " Track Momentum " << "\tMag " <<  trackMomentums[i].Mag() 
//                                       << "\tX " <<  trackMomentums[i].X()
//                                       << "\tY " <<  trackMomentums[i].Y()
//                                       << "\tZ " <<  trackMomentums[i].Z();
//   }
//   LOG(debug3) << " ====================================================== ";

//   LOG(debug3) << " =========================== Track Momentum Loss =========================== ";
//   for(size_t i = 0; i < tarckMomentumLosses.size() ; ++i)
//   {
//     LOG(debug3) << " Track Momentum " << "\tMag " <<  tarckMomentumLosses[i].Mag() 
//                                       << "\tX " <<  tarckMomentumLosses[i].X()
//                                       << "\tY " <<  tarckMomentumLosses[i].Y()
//                                       << "\tZ " <<  tarckMomentumLosses[i].Z();
//   }
//   LOG(debug3) << " ====================================================== ";

//   Bool_t isValid = (!std::isnan(momentum.Mag()) && momentum.Mag()!=inf && trackMomentums.size() > 0 );

//   return isValid;
// }


Bool_t FgdGraphGenFitRecon::CalculateInitialMomentum(const std::vector<ReconHit>& track,const TVector3& magField, TVector3& momentum , TVector3& momentumLoss)
{
  if(track.empty())
  {
      momentum.SetXYZ(0.,0.,0);
      momentumLoss.SetXYZ(0.,0.,0);
      return false;
  }
  static const Double_t inf = std::numeric_limits<Double_t>::infinity();
  TVector3 trackVector = track[track.size()-1].fHitPos - track[0].fHitPos;

  TVector3 p1 = track[0].fHitPos;
  TVector3 p2 = track[track.size()/2].fHitPos;
  TVector3 p3 = track[track.size()-1].fHitPos;

  CalculateMomentum(trackVector, p1, p2, p3, magField, momentum);

  //momentum = momentum * 1.1; // Increase momentum due to energy losses

  Bool_t isValid = (!std::isnan(momentum.Mag()) && momentum.Mag()!=inf );

  return isValid;
}


Bool_t FgdGraphGenFitRecon::CalculateMomentum(const TVector3& trackVector, const TVector3& p1, const TVector3& p2, const TVector3& p3 , const TVector3& magField, TVector3& momentum)
{
  //
  //  p [Gev/c] = e [1.6 x 10^-19 coulumb] * B [T] * R [m]
  //

  Bool_t rc(false);

  // For calculation charge is taken as 1 unit of 'e'
  const Double_t charge = 1.;

  Double_t inf = std::numeric_limits<Double_t>::infinity();

  TVector3 x_axis(1,0,0);
  TVector3 y_axis(0,1,0);
  TVector3 z_axis(0,0,1);

  Double_t x_angle = x_axis.Angle(trackVector);
  Double_t y_angle = y_axis.Angle(trackVector);
  Double_t z_angle = z_axis.Angle(trackVector);

  // For each magnetic field plane calculate it for the perpendicular projections
  if(magField.X()!=0)
  {
    TVector3 mag_point1(0,p1.Y(), p1.Z());
    TVector3 mag_point2(0,p2.Y(), p2.Z());
    TVector3 mag_point3(0,p3.Y(), p3.Z());
    Double_t radius = GetRadius(mag_point1,mag_point2,mag_point3); // radius is returned in [cm]

    if(!std::isnan(radius) && radius!=inf)
    {
      Double_t R = radius/100.; // convert in meters
      Double_t magField_T = magField.X() / 10.; // convert from kGauss to Tesla units
      Double_t mom = charge * R * magField_T;

      Double_t mom_x = std::cos(x_angle) * mom;
      Double_t mom_y = std::cos(y_angle) * mom;
      Double_t mom_z = std::cos(z_angle) * mom;

      momentum.SetX(mom_x);
      momentum.SetY(mom_y);
      momentum.SetZ(mom_z);
      rc = true;
    }
  }

  if(magField.Y()!=0)
  {
    TVector3 mag_point1(p1.X(), 0. , p1.Z());
    TVector3 mag_point2(p2.X(), 0. , p2.Z());
    TVector3 mag_point3(p3.X(), 0. , p3.Z());
    Double_t radius = GetRadius(mag_point1,mag_point2,mag_point3); // radius is returned in [cm]

    if(!std::isnan(radius) && radius!=inf)
    {
      Double_t R = radius/100.; // convert in meters
      Double_t magField_T = magField.Y() / 10.; // convert from kGauss to Tesla units
      Double_t mom = charge * R * magField_T;

      Double_t mom_x = std::cos(x_angle) * mom + momentum.X();
      Double_t mom_y = std::cos(y_angle) * mom + momentum.Y();
      Double_t mom_z = std::cos(z_angle) * mom + momentum.Z();

      momentum.SetX(mom_x);
      momentum.SetY(mom_y);
      momentum.SetZ(mom_z);
      rc = true;
    }
  }
  
  if(magField.Z()!=0)
  {
    TVector3 mag_point1(p1.X(), p1.Y() , 0.);
    TVector3 mag_point2(p2.X(), p2.Y() , 0.);
    TVector3 mag_point3(p3.X(), p3.Y() , 0.);
    Double_t radius = GetRadius(mag_point1,mag_point2,mag_point3); // radius is returned in [cm]

    if(!std::isnan(radius) && radius!=inf)
    {
      Double_t R = radius/100.; // convert in meters
      Double_t magField_T = magField.Z() / 10.; // convert from kGauss to Tesla units
      Double_t mom = charge * R * magField_T;

      Double_t mom_x = std::cos(x_angle) * mom + momentum.X();
      Double_t mom_y = std::cos(y_angle) * mom + momentum.Y();
      Double_t mom_z = std::cos(z_angle) * mom + momentum.Z();

      momentum.SetX(mom_x);
      momentum.SetY(mom_y);
      momentum.SetZ(mom_z);
      rc = true;
    }
  }

  return rc;
}




Bool_t FgdGraphGenFitRecon::CalculateCalorimetricMomentum(const std::vector<ReconHit>& track, TVector3& momentum)
{
  static const double EdepToPhotConv_FGD = 70.8 / CLHEP::MeV; // contains both collection in fiber and edep->gamma conversion 
  static const double EdepToPhotConv_SuperFGD = EdepToPhotConv_FGD * 1.3;

  static const Double_t a = -22.67;
  static const Double_t b =  12.124;

  static const TVector3 x_axis(1,0,0);
  static const TVector3 y_axis(0,1,0);
  static const TVector3 z_axis(0,0,1);

  if(track.empty()) return false;

  Double_t totalEdep(0.);

  for(size_t i = 0; i<track.size(); ++i)
  {
    const ReconHit& hit = track[i];
    totalEdep+= hit.fphotons.X() + hit.fphotons.Y() +hit.fphotons.Z();
  }

  Double_t trackLength = track.size() * f_step_X;

  LOG(info) << "Cal totalEdep-> " << totalEdep;
  totalEdep = 10*totalEdep/EdepToPhotConv_FGD;
  LOG(info) << "Cal totalEdep [after]-> " << totalEdep;
  Double_t dedx1 = totalEdep/trackLength;

  Double_t dedx = RevertToDeDx(track);

  Double_t mom = (dedx - b)/a;

  LOG(info) << "Cal dedx-> " << dedx;
  LOG(info) << "Cal dedx1-> " << dedx1;
  LOG(info) << "Cal dedx with coeff -> " << dedxToP(dedx);
  LOG(info) << "Cal trackLength-> " << trackLength;

  TVector3 trackVector = track[track.size()-1].fHitPos - track[0].fHitPos;
  Double_t x_angle = x_axis.Angle(trackVector);
  Double_t y_angle = y_axis.Angle(trackVector);
  Double_t z_angle = z_axis.Angle(trackVector);


  Double_t mom_x = std::cos(x_angle) * mom;
  Double_t mom_y = std::cos(y_angle) * mom;
  Double_t mom_z = std::cos(z_angle) * mom;

  momentum.SetX(mom_x);
  momentum.SetY(mom_y);
  momentum.SetZ(mom_z);
  
  return true;
}

// Calculate the radius from 3 points using the "Menger curvature" theorem
Double_t FgdGraphGenFitRecon::GetRadius(const TVector3& p1, const TVector3& p2, const TVector3& p3)
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

Int_t FgdGraphGenFitRecon::GetPdgCode(const TVector3& momentum, const TVector3& momentumLoss , Int_t tryFit, Int_t& momCoeff)
{
  //Int_t pdgCode = 13; // Muon pdg code

  // for(size_t k=0; k<fpdgFromMomLoss.size(); ++k)
  // {
  //   PdgFromMomentumLoss& p = fpdgFromMomLoss[k];
  //   Double_t momtoLoss = momentum.Mag()/ momentumLoss.Mag();
  //   if(p.GetPdg(momtoLoss, pdgCode))
  //   {
  //     break;
  //   }
  // }

  Int_t pdgGuess = tryFit / 3;

  Int_t pdgCode = PDG_FIT_TRACK::MUON;
  switch(pdgGuess)
  {
      case 0: 
          {
            pdgCode = PDG_FIT_TRACK::MUON;
          }
        break;

      case 1:
          {
            pdgCode = PDG_FIT_TRACK::PION;
          }
        break;

      case 2:
          {
            pdgCode = PDG_FIT_TRACK::PROTON;
          }
        break;

      case 3:
          {
            pdgCode = PDG_FIT_TRACK::ELECTRON;
          }
        break;

      default:
          {

          }
        break;
  }

  momCoeff = (1 << (tryFit % 3));
  //LOG(debug) << "Trying with pdg code = " << pdgCode << " , momCoeff = " << momCoeff;
  return pdgCode;
}

void FgdGraphGenFitRecon::ConvertHitToVec(std::vector<TVector3>& points, std::vector<ReconHit>& hits)
{
  for(size_t i=0; i< hits.size(); ++i)
  {
    points.emplace_back(hits[i].fHitPos.X()
                      , hits[i].fHitPos.Y() 
                      , hits[i].fHitPos.Z());
  }
}

void FgdGraphGenFitRecon::FitTracks(std::vector<std::vector<ReconHit>>& foundTracks)
{
    fTracksArray->Delete();
    
    // init fitter
    std::shared_ptr<genfit::AbsKalmanFitter> fitter = make_shared<genfit::KalmanFitterRefTrack>();
    //std::shared_ptr<genfit::AbsKalmanFitter> fitter = make_shared<genfit::KalmanFitter>();
    fitter->setMinIterations(fminGenFitInterations);
    fitter->setMaxIterations(fmaxGenFitIterations);
    fitter->setDebugLvl(fDebuglvl_genfit);

    std::vector<ReconHit>& longestTrack = foundTracks[0];
    size_t ltId(0);
    for(size_t i = 1; i <  foundTracks.size() ; ++i)
    {
      if(foundTracks[i].size() > longestTrack.size())
      {
          longestTrack = foundTracks[i];
          ltId = i;
      }
    }

    FgdTMVAEventRecord& tvmaEventRecord = feventRecords[feventNum];
    // Assume that the longest track is that of the muon
    TVector3 mom(0,0,0);
    TVector3 calMom(0,0,0);
    if(FitTrack(longestTrack,fitter, PDG_FIT_TRACK::MUON, mom, calMom, ltId))
    {
      tvmaEventRecord.SetMuonFitMom(mom);
    }
    tvmaEventRecord.SetMuonCalorimetricMom(calMom);
}

bool FgdGraphGenFitRecon::FitTrack(std::vector<ReconHit>& track
    , std::shared_ptr<genfit::AbsKalmanFitter>& fitter
    , Int_t pdg
    , TVector3& momentum
    , TVector3& calMom
    , int trackId)
{
    TVector3 magField = fgdConstructor.GetMagneticField(); // values are in kGauss
    std::vector<genfit::Track*> genTracks;
    static int detId(1); // Detector id, it is the same, we only have one detector
    Bool_t rc(false);

    std::vector<ReconHit>& hitsOnTrack = track;

    // Set lower limit on track size
    if(hitsOnTrack.size()<fminHits)
    {
      LOG(debug) << "Track " << trackId << " below limit, continue with next track (" << hitsOnTrack.size() << " < " << fminHits << ")";
      return false;
    }
    
  
    TVector3 posM(hitsOnTrack[0].fHitPos.X(),hitsOnTrack[0].fHitPos.Y(),hitsOnTrack[0].fHitPos.Z());
    TVector3 momM(0,0,0);
    //TVector3 momLoss(0,0,0);

    //momM = hitsOnTrack[0].fmom;

    // if(!CalculateInitialMomentum(hitsOnTrack, magField, momM, momLoss))
    // {
    //   //LOG(debug) << "Track " << trackId << " unable to extract momentum. Continue with next track";
    //   //return false;
    // }

    calMom = getCalorimetricMomentum(hitsOnTrack);
    momM = calMom;

    /// TODO: Test template fit
    // LOG(info) << "Before fit ";
    // TVector3 momtum;
    // bool res = freconTemplate.FitTrack(track
    //                     ,  true
    //                     ,  pdg
    //                     ,  posM
    //                     ,  momM
    //                     ,  false
    //                     ,  momtum);

    // LOG(info) << "Fitted Momentum " << "[" << (res?"true":"false") <<"] " << momtum.Mag();

    // if(CalculateCalorimetricMomentum(hitsOnTrack, calMom))
    // {
    //     LOG(debug) << " \tCalorimetricMomentum [" << calMom.Mag() << "]" << "(" << calMom.X() << "," << calMom.Y() << "," << calMom.Z() << ")";
    // }

    // approximate covariance
    const double resolution = 2;// Default in example is 0.1;
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
    //stateSmeared.setPosMomCov(posM, calMom, covM);
    
    // create track
    TVectorD seedState(6);
    TMatrixDSym seedCov(6);
    stateSmeared.get6DStateCov(seedState, seedCov);

    genfit::Track* toFitTrack = new genfit::Track(rep, seedState, seedCov);
    // genfit::Track* toFitTrack = new genfit::Track();
    // toFitTrack->setStateSeed(seedState);
    // toFitTrack->setCovSeed(seedCov);
    // toFitTrack->addTrackRep(rep);

    size_t bhlimit = hitsOnTrack.size();

    LOG(debug) << "******************************************* ";
    LOG(debug) << "******    Track "<< trackId << "  ************************";
    LOG(debug) << "******************************************* ";
    LOG(debug) << " \tPdg code " << pdg << " [MC pdg = " << hitsOnTrack[0].fpdg <<" ]";
    LOG(debug) << " \tHits in track "<< bhlimit;
    LOG(debug) << " \tMC Momentum [" << hitsOnTrack[0].fmom.Mag() << "]" << "(" << hitsOnTrack[0].fmom.X() << "," << hitsOnTrack[0].fmom.Y() << "," << hitsOnTrack[0].fmom.Z() << ")";
    LOG(debug) << " \tEstimated Momentum [" << momM.Mag() << "]" << "(" << momM.X() << "," << momM.Y() << "," << momM.Z() << ")";
    //LOG(debug) << " \tTrack Momentum Loss [" << momLoss.Mag() << "]" << "(" << momLoss.X() << "," << momLoss.Y() << "," << momLoss.Z() << ")";
    //LOG(debug) << " \tMomentum / Momentum Loss [" << momM.Mag()/momLoss.Mag() << "]";
    LOG(debug) << " \tCalorimetric Momentum [" << calMom.Mag() << "]" << "(" << calMom.X() << "," << calMom.Y() << "," << calMom.Z() << ")";
    
    FgdTMVAEventRecord& tvmaEventRecord = feventRecords[feventNum];
    const std::vector<std::pair<Int_t, TVector3>>& particles = tvmaEventRecord.GetPrimaryParticles();
    for(size_t p = 0; p < particles.size(); ++p)
    {
      std::pair<Int_t, TVector3> pp = particles[p];
      if(genie::pdg::IsMuon(pp.first) || genie::pdg::IsAntiMuon(pp.first))
      {
        LOG(debug) << " \tEvent Muon Momentum [" << pp.second.Mag() << "]" << "(" << pp.second.X() << "," << pp.second.Y() << "," << pp.second.Z() << ")";
        break;
      }
    }
    
    std::vector<genfit::AbsMeasurement*> measurements;
    for(Int_t bh = 0; bh < bhlimit; ++bh)
    //for(Int_t bh = 0; bh < bhlimit; bh+=2)
    {
      TVectorD hitPos(3);
      hitPos(0) = fuseSmoothPos? hitsOnTrack[bh].fsmoothco.X() : hitsOnTrack[bh].fHitPos.X();
      hitPos(1) = fuseSmoothPos? hitsOnTrack[bh].fsmoothco.Y() : hitsOnTrack[bh].fHitPos.Y();
      hitPos(2) = fuseSmoothPos? hitsOnTrack[bh].fsmoothco.Z() : hitsOnTrack[bh].fHitPos.Z();

      genfit::AbsMeasurement* measurement = new genfit::SpacepointMeasurement(hitPos, hitCov, detId, bh, nullptr);
      measurements.emplace_back(measurement);
    }
    toFitTrack->insertPoint(new genfit::TrackPoint(measurements, toFitTrack));

    try
    {
      //((genfit::KalmanFitter*)fitter.get())->useSquareRootFormalism();

      // int nFailedHits(0);
      // bool res = ((genfit::KalmanFitterRefTrack*)fitter.get())->prepareTrack(toFitTrack, rep, true, nFailedHits);
      // LOG(debug) <<"prepareTrack " << res;
      // while(res)
      // {
      //   res = ((genfit::KalmanFitterRefTrack*)fitter.get())->prepareTrack(toFitTrack, rep, true, nFailedHits);
      //   LOG(debug) <<"prepareTrack " << res;
      // }
    
      //check
      toFitTrack->checkConsistency();

      // do the fit
      fitter->processTrack(toFitTrack, true);
      //fitter->processTrackWithRep(toFitTrack, toFitTrack->getCardinalRep(), false);
      // double chi2(0);
      // double ndf(0);
      // int direction(-1);
      // ((genfit::KalmanFitterRefTrack*)fitter.get())->fitTrack(toFitTrack, rep, chi2, ndf, direction);

      //check
      toFitTrack->checkConsistency();

      PrintFitTrack(*toFitTrack);

      LOG(debug) <<"******************************************* ";
      genfit::FitStatus* fiStatuStatus = toFitTrack->getFitStatus(rep);
      rc = fiStatuStatus->isFitted() && (fiStatuStatus->isFitConverged() || fiStatuStatus->isFitConvergedPartially());

      if(rc)
      {
        genTracks.push_back(toFitTrack);
      }
      
      const genfit::MeasuredStateOnPlane& me = (*toFitTrack).getFittedState();
      momentum = me.getMom();
      LOG(debug) << " \tFitted Momentum [" << momentum.Mag() << "]" << "(" << momentum.X() << "," << momentum.Y() << "," << momentum.Z() << ")";
    }
    catch(genfit::Exception& e)
    {
        LOG(error) <<"Exception, when tryng to fit track";
        LOG(error) << "e.what() " << e.what();
        LOG(error) << "e.getExcString() " << e.getExcString();
    }
  
    if(isGenFitVisualization)
    {
      fdisplay->addEvent(genTracks);
    }

    return rc;
}

void FgdGraphGenFitRecon::ExtractTMVAdata(std::vector<ReconHit>& allhits)
{
    if(feventNum >= feventRecords.size())
    {
        LOG(fatal) << "EventData reconrds are less than the passed events!";
        throw "EventData reconrds are less than the passed events!";
    }

    FgdTMVAEventRecord& tvmaEventRecord = feventRecords[feventNum];

    Int_t sumTotalCubes = 0;
    TVector3 sumTotalPhoto(0,0,0);

    for(size_t i = 0; i <  allhits.size() ; ++i)
    {
        ReconHit& hit = allhits[i];
        sumTotalPhoto +=hit.fphotons; 
    }

    tvmaEventRecord.SetTotalPhotons(sumTotalPhoto);
    tvmaEventRecord.SetTotalCubes(allhits.size());
}

void FgdGraphGenFitRecon::DefineMaterials() 
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

void FgdGraphGenFitRecon::PrintFitTrack(genfit::Track& fitTrack)
{
  const genfit::MeasuredStateOnPlane& me = fitTrack.getFittedState();
  //LOG(debug)<< "Momentum  " << (me.getMom()).Mag();
  //LOG(debug)<< " X  " << (me.getMom()).X()<< " Y " << (me.getMom()).Y()<< " Z  " << (me.getMom()).Z();
  LOG(debug)<< "\tFitted Momentum  [" << (me.getMom()).Mag() <<"]" << "(" << (me.getMom()).X() << "," << (me.getMom()).Y() << "," << (me.getMom()).Z() << ")";
  genfit::FitStatus* fiStatuStatus = fitTrack.getFitStatus();
  fiStatuStatus->Print();

  LOG(debug)<< "fiStatuStatus->isFitted()  " << fiStatuStatus->isFitted();
  LOG(debug)<< "fiStatuStatus->isFitConverged()  " << fiStatuStatus->isFitConverged();
  LOG(debug)<< "fiStatuStatus->isFitConvergedFully()  " << fiStatuStatus->isFitConvergedFully();
  LOG(debug)<< "fiStatuStatus->isFitConvergedPartially()  " << fiStatuStatus->isFitConvergedPartially();
  LOG(debug)<< "fitTrack.getNumPoints() " << fitTrack.getNumPoints();
  LOG(debug)<< "fitTrack.getNumReps() " << fitTrack.getNumReps();
  
}

Long_t FgdGraphGenFitRecon::ArrInd(int x, int y, int z)
{
  return (x*f_bin_Y*f_bin_Z + y*f_bin_Z+z);
}

// -------------------------------------------------------------------------

}// namespace superfgd
}// namespace reconstruction
}// namespace esbroot
