#include "EsbReconstruction/EsbSuperFGD/FgdTMVAData.h"
#include "EsbReconstruction/EsbSuperFGD/FgdReconTemplate.h"
//#include "EsbData/EsbSuperFGD/FgdDetectorPoint.h"
#include "EsbDigitizer/EsbSuperFGD/FgdDigitizer.h"

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
#include <TTree.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>

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

#define MAX_LENGTH_TRACKS_TO_RECORD 3 // Number of cubes

// -----   Default constructor   -------------------------------------------
FgdTMVAData::FgdTMVAData() : FgdMCGenFitRecon(), feventNum(0), fmagField_X(0.), fmagField_Y(0.), fmagField_Z(0.)
    , fMaxtrack(1),fMaxTotph(1),fMaxCubes(1), fMaxTrph(1) , fMaxnuE(-1), fMaxTrueEdep(-1)
{
}
// -------------------------------------------------------------------------

// -----   Constructor   -------------------------------------------
FgdTMVAData::FgdTMVAData(const char* name
                          , const char* geoConfigFile
                          , const char* mediaFile
                          , const char* eventData
                          , const char* outputRootFile
                          , Int_t verbose
                          , double debugLlv) :
  FgdMCGenFitRecon(name, geoConfigFile, mediaFile, verbose, 
                    debugLlv, false /* no visualization */, "D")
    , feventData(eventData), foutputRootFile(outputRootFile)
    , feventNum(0), fmagField_X(0.), fmagField_Y(0.), fmagField_Z(0.)
    , fMaxtrack(1),fMaxTotph(1),fMaxCubes(1), fMaxTrph(1), fMaxTotEdep(1)
    , fMaxTotPe(1)
{ 
    fpdgDB = make_shared<TDatabasePDG>();
}
// -------------------------------------------------------------------------



// -----   Destructor   ----------------------------------------------------
FgdTMVAData::~FgdTMVAData() 
{
    if(fTracksArray) 
    {
        fTracksArray->Delete();
        delete fTracksArray;
    } 
}
// -------------------------------------------------------------------------



// -----   Public method Init   --------------------------------------------
InitStatus FgdTMVAData::Init() 
{   
    FgdMCGenFitRecon::Init();

    fmagField_X = fParams.ParamAsDouble(esbroot::geometry::superfgd::DP::magField_X);
    fmagField_Y = fParams.ParamAsDouble(esbroot::geometry::superfgd::DP::magField_Y);
    fmagField_Z = fParams.ParamAsDouble(esbroot::geometry::superfgd::DP::magField_Z);

    std::ifstream eventFileStream;
    try
    {        
        eventFileStream.open(feventData.c_str(), std::ios::in);
        Int_t id = 0;

        if(eventFileStream.is_open())
        {
            std::string line;
            while(std::getline(eventFileStream,line))
            {
                feventRecords.emplace_back(FgdTMVAEventRecord(line));
                fhitCoordinates[id] = std::vector<TVector3>{};
                fhitPhotons[id] = std::vector<TVector3>{};
                ++id;
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

void FgdTMVAData::OutputFileInit(FairRootManager* manager)
{
    FgdMCGenFitRecon::OutputFileInit(manager);
}

// -------------------------------------------------------------------------



// -----   Public methods   --------------------------------------------

void FgdTMVAData::Exec(Option_t* opt) 
{  
  try
  {
    std::vector<ReconHit> allhits;
    std::vector<std::vector<ReconHit>> foundTracks;

    bool rc = GetHits(allhits);

    if(rc)
    { 
      LOG(debug) <<" Hits to retrieve stats from " << allhits.size();
      SplitTrack(allhits, foundTracks);
    }

    if(rc)
    {
      LOG(debug) <<" Found tracks to process " << foundTracks.size();
      ProcessStats(foundTracks);
    }
    else
    {
      LOG(error) << " Could not find clean hits or tracks! ";
      ++feventNum;
    }
  }
  catch(genfit::Exception& e)
  {
      LOG(fatal) << "Exception, when tryng to fit track";
      LOG(fatal) << e.what();
  }
}

Bool_t FgdTMVAData::ProcessStats(std::vector<std::vector<ReconHit>>& foundTracks)
{
    if(feventNum >= feventRecords.size())
    {
        LOG(fatal) << "EventData reconrds are less than the passed events!";
        throw "EventData reconrds are less than the passed events!";
    }

    FgdTMVAEventRecord& tvmaEventRecord = feventRecords[feventNum];
    tvmaEventRecord.ReadEventData();

    if(fMaxnuE < tvmaEventRecord.GetNuE())
    {
        fMaxnuE = tvmaEventRecord.GetNuE();
    }

    for(size_t i = 0; i <  foundTracks.size() ; ++i)
    {
        std::vector<ReconHit>& hitsOnTrack = foundTracks[i];
        // Sort by time, the 1st hit in time is the start of the track
        std::sort(hitsOnTrack.begin(), hitsOnTrack.end(), [](ReconHit& bh1, ReconHit& bh2){return bh1.ftime<bh2.ftime;});
    }


    // 1. Extract data for all hits and total cubes hits, total photons
    Int_t sumTotalCubes = 0;
    TVector3 sumTotalPhoto(0,0,0);
    Float_t sumEdep = 0;
    Float_t sumPe = 0;
    Float_t sumTrueEdep = 0;
    for(size_t i = 0; i <  foundTracks.size() ; ++i)
    {
        std::vector<ReconHit>& hitsOnTrack = foundTracks[i];
        if(hitsOnTrack.empty() || !isParticleAllowed(hitsOnTrack[0].fpdg)) continue;

        for(size_t j = 0; j < hitsOnTrack.size(); ++j)
        {   
            ReconHit& hit = hitsOnTrack[j];
            sumTotalPhoto +=hit.fphotons;
            sumTotalCubes++;
            sumEdep += CalculatePhotoEdep(hit);
            sumPe += hit.fpe;
            sumTrueEdep += hit.fEdep;

            fhitCoordinates[feventNum].emplace_back(hit.fmppcLoc);
            fhitPhotons[feventNum].emplace_back(hit.fphotons);
        }   
    }

    tvmaEventRecord.SetTotalPhotons(sumTotalPhoto);
    tvmaEventRecord.SetTotalCubes(sumTotalCubes); 
    tvmaEventRecord.SetTotalEdep(sumEdep);
    tvmaEventRecord.SetTrueEdep(sumTrueEdep);
    tvmaEventRecord.SetPe(sumPe);
    Float_t sumTotph = sumTotalPhoto.X() + sumTotalPhoto.Y() + sumTotalPhoto.Z();
    if(fMaxTotph < sumTotph)
    {
        fMaxTotph = sumTotph;
    }
    if(fMaxCubes < sumTotalCubes)
    {
        fMaxCubes = sumTotalCubes;
    }
    if(fMaxTotEdep < sumEdep)
    {
        fMaxTotEdep = sumEdep;
    }
    if(fMaxTotPe < sumPe)
    {
        fMaxTotPe = sumPe;
    }
    if(fMaxTrueEdep < sumTrueEdep)
    {
        fMaxTrueEdep = sumTrueEdep;
    }

    


    // 2. Get data only for the first three longest tracks
    Int_t min_Track_Length = fParams.ParamAsInt(esbroot::geometry::superfgd::DP::FGD_MIN_TRACK_LENGTH);
    std::sort(foundTracks.begin(), foundTracks.end(), 
        [](std::vector<ReconHit>& tr1, std::vector<ReconHit>& tr2){return tr1.size()>tr2.size();});

    size_t limit = foundTracks.size() > MAX_LENGTH_TRACKS_TO_RECORD? MAX_LENGTH_TRACKS_TO_RECORD: foundTracks.size();
    std::vector<Float_t> track = {0.,0.,0.};
    std::vector<Float_t> phototrack= {0.,0.,0.};
    for(size_t i = 0; i <  limit ; ++i)
    {
        std::vector<ReconHit>& hitsOnTrack = foundTracks[i];
        if(hitsOnTrack.size() < min_Track_Length) break;

        // size_t trackLenght = hitsOnTrack.size();
        
        // if(fMaxtrack < trackLenght)
        // {
        //     fMaxtrack = trackLenght;
        // }
        // track[i] = trackLenght;

        Float_t trackLenght = hitsOnTrack[hitsOnTrack.size()-1].ftrackLengthOrigin;
        if(fMaxtrack < trackLenght)
        {
            fMaxtrack = trackLenght;
        }
        track[i] = trackLenght;

        Float_t sumPh(0.);
        for(size_t j = 0; j < hitsOnTrack.size(); ++j)
        {   
            ReconHit& hit = hitsOnTrack[j];
            sumPh +=hit.fphotons.X() + hit.fphotons.Y() + hit.fphotons.Z();
        }
        if(fMaxTrph < sumPh)
        {
            fMaxTrph = sumPh;
        }
        phototrack[i] = sumPh;

    }
    ftrackLenghts.emplace_back(std::move(track));
    ftrackPhotos.emplace_back(std::move(phototrack));

    feventRecords[feventNum].SetHasHits(true);
    ++feventNum; // Increment to next event from eventData read from simulation`s genie export
}

Double_t FgdTMVAData::CalculatePhotoEdep(ReconHit& hit)
{
    using namespace esbroot::digitizer::superfgd;

    TVector3& pe_1_dir = hit.fph1;
    TVector3& dist_1 = hit.fmppc1;
    TVector3& pe_2_dir = hit.fph2;
    TVector3& dist_2 = hit.fmppc2;

    double& time = hit.ftime;
    double charge = 1.0; // not used, but is a parameter from legacy


    double pe_1_x = FgdDigitizer::RevertyMPPCResponse(pe_1_dir.X());
    double pe_1_y = FgdDigitizer::RevertyMPPCResponse(pe_1_dir.Y());
    double pe_1_Z = FgdDigitizer::RevertyMPPCResponse(pe_1_dir.Z());

    double pe_2_x = FgdDigitizer::RevertyMPPCResponse(pe_2_dir.X());
    double pe_2_y = FgdDigitizer::RevertyMPPCResponse(pe_2_dir.Y());
    double pe_2_z = FgdDigitizer::RevertyMPPCResponse(pe_2_dir.Z());

    FgdDigitizer::RevertFiberResponse(pe_1_x, time, dist_1.X());
    FgdDigitizer::RevertFiberResponse(pe_1_y, time, dist_1.Y());
    FgdDigitizer::RevertFiberResponse(pe_1_Z, time, dist_1.Z());

    FgdDigitizer::RevertFiberResponse(pe_2_x, time, dist_2.X());
    FgdDigitizer::RevertFiberResponse(pe_2_y, time, dist_2.Y());
    FgdDigitizer::RevertFiberResponse(pe_2_z, time, dist_2.Z());

    // Deposited energy hit.fEdep and tracklength hit.ftrackLength
    // are used to calculate the CBIRKS coefficients together with dedx (energy losses)
    // to be able to revert from the photons

    Double_t x_1 = FgdDigitizer::RevertScintiResponse(hit.fEdep, hit.ftrackLength, charge, pe_1_x);
    Double_t y_1 = FgdDigitizer::RevertScintiResponse(hit.fEdep, hit.ftrackLength, charge, pe_1_y);
    Double_t z_1 = FgdDigitizer::RevertScintiResponse(hit.fEdep, hit.ftrackLength, charge, pe_1_Z);

    Double_t x_2 = FgdDigitizer::RevertScintiResponse(hit.fEdep, hit.ftrackLength, charge, pe_2_x);
    Double_t y_2 = FgdDigitizer::RevertScintiResponse(hit.fEdep, hit.ftrackLength, charge, pe_2_y);
    Double_t z_2 = FgdDigitizer::RevertScintiResponse(hit.fEdep, hit.ftrackLength, charge, pe_2_z);

    Double_t x = x_1 + x_2;
    Double_t y = y_1 + y_2;
    Double_t z = z_1 + z_2;

    Double_t totalEdep =  x + y + z;
    return totalEdep;
}


void FgdTMVAData::FinishTask()
{
    TFile * outFile = new TFile(foutputRootFile.c_str(), "RECREATE", "TVMA data from Fgd Detector");
    outFile->SetCompressionLevel(9);

    // 1. Write event data, hit positions and total photons of hits
    // FgdTMVAEventRecord* data = nullptr;
    // TClonesArray* hitCoordinates = new TClonesArray(TVector3::Class());
    // TClonesArray& hitcref = *hitCoordinates;

    // TClonesArray* hitPhotons = new TClonesArray(TVector3::Class());
    // TClonesArray& photoref = *hitPhotons;

    // TTree * outTree = new TTree(esbroot::geometry::superfgd::DP::FGD_TMVA_DATA_TTREE.c_str()
    //                             ,esbroot::geometry::superfgd::DP::FGD_TMVA_DATA_ROOT_FILE.c_str());
    // outTree->Branch(esbroot::geometry::superfgd::DP::FGD_TMVA_DATA_BRANCH.c_str(), &data);
    // outTree->Branch(esbroot::geometry::superfgd::DP::FGD_TMVA_HIT_ARRAY_BRANCH.c_str(), &hitCoordinates);
    // outTree->Branch(esbroot::geometry::superfgd::DP::FGD_TMVA_PHOTO_ARRAY_BRANCH.c_str(), &hitPhotons);

    // for(size_t ind = 0 ; ind < feventRecords.size(); ind++)
    // {
    //     data = &feventRecords[ind];

    //     // 1. Copy hitposition
    //     for(Int_t i = 0 ; i < fhitCoordinates[ind].size(); i++)
    //     {
    //         new(hitcref[i]) TVector3(fhitCoordinates[ind][i]);
    //     }

    //     // 2. Copy photons
    //     for(Int_t i = 0 ; i < fhitPhotons[ind].size(); i++)
    //     {
    //         new(photoref[i]) TVector3(fhitPhotons[ind][i]);
    //     }

    //     outTree->Fill();

    //     hitCoordinates->Clear();
    // }
    // outFile->WriteTObject(outTree);
    // =================================================================

    // 2. Write simple format for analysis
    // Containing total photons and nu energy
    // Float_t totalPhX = 0;
    // Float_t totalPhY = 0;
    // Float_t totalPhZ = 0;
    // Float_t totalCubes = 0;
    // Float_t nuE = 0.;
    // Bool_t isCC = false;
    // Bool_t isQuasiE = false;
    // Float_t nuPdgph = 0;
    // TTree * totalPhTree = new TTree(esbroot::geometry::superfgd::DP::FGD_TOTAL_PHOTONS_TTREE.c_str()
    //                             ,esbroot::geometry::superfgd::DP::FGD_TMVA_DATA_ROOT_FILE.c_str());
    // totalPhTree->Branch("totalPhotonsX", &totalPhX);
    // totalPhTree->Branch("totalPhotonsY", &totalPhY);
    // totalPhTree->Branch("totalPhotonsZ", &totalPhZ);
    // totalPhTree->Branch("totalCubes", &totalCubes);
    // totalPhTree->Branch("nuEnergy", &nuE);
    // totalPhTree->Branch("isCC", &isCC);
    // totalPhTree->Branch("isQuasiE", &isQuasiE);
    // totalPhTree->Branch("nuPdg", &nuPdgph);
    // totalPhTree->Branch("magFieldX", &fmagField_X);
    // totalPhTree->Branch("magFieldY", &fmagField_Y);
    // totalPhTree->Branch("magFieldZ", &fmagField_Z);
    
    // for(size_t ind = 0 ; ind < feventRecords.size(); ind++)
    // {
    //     data = &feventRecords[ind];
    //     totalPhX = data->GetTotalPhotons().X();
    //     totalPhY = data->GetTotalPhotons().Y();
    //     totalPhZ = data->GetTotalPhotons().Z();
    //     totalCubes = data->GetTotalCubes();
    //     nuE = data->GetNuE();
    //     isCC = data->IsWeakCC();
    //     isQuasiE = data->IsQuasiElastic();
    //     nuPdgph = data->GetNuPdg();

    //     totalPhTree->Fill();
    // }
    // outFile->WriteTObject(totalPhTree);
    // =================================================================


    // 3. Write track projections
    // Containing total photons and nu energy
    // TTree * trackPrjTree = new TTree(esbroot::geometry::superfgd::DP::FGD_TRACK_PROJECTION_TTREE.c_str()
    //                             ,esbroot::geometry::superfgd::DP::FGD_TMVA_DATA_ROOT_FILE.c_str());

    // std::vector<Int_t> x_projections(f_bin_X+1,0);
    // std::vector<Int_t> y_projections(f_bin_Y+1,0);
    // std::vector<Int_t> z_projections(f_bin_Z+1,0);
    // Float_t nuEnergy = 0.;
    // Float_t nuPdg = 0;
    // std::stringstream ss;

    // for(size_t x = 0; x < f_bin_X; ++x)
    // {
    //     ss << "x"<< x;
    //     trackPrjTree->Branch(ss.str().c_str(), &x_projections[x]);
    //     ss.str("");
    // }

    // for(size_t y = 0; y < f_bin_Y; ++y)
    // {
    //     ss << "y"<< y;
    //     trackPrjTree->Branch(ss.str().c_str(), &y_projections[y]);
    //     ss.str("");
    // }

    // for(size_t z = 0; z < f_bin_Z; ++z)
    // {
    //     ss << "z"<< z;
    //     trackPrjTree->Branch(ss.str().c_str(), &z_projections[z]);
    //     ss.str("");
    // }

    // trackPrjTree->Branch("nuEnergy", &nuEnergy);
    // trackPrjTree->Branch("nuPdg", &nuPdg);

    // const Int_t events = feventRecords.size();
    // for(size_t ind = 0 ; ind < events; ind++)
    // {
    //     data = &feventRecords[ind];
    //     nuEnergy = data->GetNuE();
    //     nuPdg = data->GetNuPdg();
    //     const Int_t cubes = fhitCoordinates[ind].size();
    //     for(Int_t i = 0 ; i < cubes; i++)
    //     {
    //         TVector3& coordinate = fhitCoordinates[ind][i];
    //         TVector3& photons = fhitPhotons[ind][i];

    //         x_projections[(Int_t)coordinate.X()] = photons.X();
    //         y_projections[(Int_t)coordinate.Y()] = photons.Y();
    //         z_projections[(Int_t)coordinate.Z()] = photons.Z();
    //     }
    //     trackPrjTree->Fill();

    //     for(size_t x = 0; x < f_bin_X; ++x) x_projections[x]=0;
    //     for(size_t y = 0; y < f_bin_Y; ++y) y_projections[y]=0;
    //     for(size_t z = 0; z < f_bin_Z; ++z) z_projections[z]=0;
    // }
    // outFile->WriteTObject(trackPrjTree);                      
    // =================================================================




    // 4. Write 3 longest track projections
    // Containing total photons and nu energy
    TTree * longestTrackPrjTree = new TTree(esbroot::geometry::superfgd::DP::FGD_LONGEST_TRACK_PROJECTION_TTREE.c_str()
                                ,esbroot::geometry::superfgd::DP::FGD_TMVA_DATA_ROOT_FILE.c_str());

    longestTrackPrjTree->Branch("magFieldX", &fmagField_X);
    longestTrackPrjTree->Branch("magFieldY", &fmagField_Y);
    longestTrackPrjTree->Branch("magFieldZ", &fmagField_Z);

    Float_t tr1(0);
    Float_t tr2(0);
    Float_t tr3(0);

    longestTrackPrjTree->Branch("tr1", &tr1);
    longestTrackPrjTree->Branch("tr2", &tr2);
    longestTrackPrjTree->Branch("tr3", &tr3);

    std::vector<Float_t> ph_trVec= {0,0,0};

    Float_t ph_tr1(0);
    Float_t ph_tr2(0);
    Float_t ph_tr3(0);

    longestTrackPrjTree->Branch("ph_tr1", &ph_tr1);
    longestTrackPrjTree->Branch("ph_tr2", &ph_tr2);
    longestTrackPrjTree->Branch("ph_tr3", &ph_tr3);

    Float_t lnuEnergy = 0.;
    Float_t lnuPdg = 0;
    Float_t totCubes = 0;
    Float_t totPh = 0.;
    Float_t totalEdep = 0.;
    Float_t trueEdep = 0.;
    Float_t totalPe = 0.;

    longestTrackPrjTree->Branch("totalCubes", &totCubes);
    longestTrackPrjTree->Branch("totalPhotons", &totPh);
    longestTrackPrjTree->Branch("totalEdep", &totalEdep);
    longestTrackPrjTree->Branch("trueEdep", &trueEdep);
    longestTrackPrjTree->Branch("totalPe", &totalPe);
    longestTrackPrjTree->Branch("nuEnergy", &lnuEnergy);
    longestTrackPrjTree->Branch("nuPdg", &lnuPdg);

    const Int_t evInd = feventRecords.size();
    FgdTMVAEventRecord* dataEvent = nullptr;
    for(size_t ind = 0 ; ind < evInd && ind < feventNum; ind++)
    {
        dataEvent = &feventRecords[ind];

        bool isQuasiCC = dataEvent->IsWeakCC() && dataEvent->IsQuasiElastic();
        if(!isQuasiCC)
        {
            continue;
        }

        lnuEnergy = dataEvent->GetNuE();
        lnuPdg = dataEvent->GetNuPdg();
        totalEdep = dataEvent->GetTotalEdep();
        totalPe = dataEvent->GetPe();
        trueEdep = dataEvent->GetTrueEdep();

        totPh = dataEvent->GetTotalPhotons().X() + dataEvent->GetTotalPhotons().Y() + dataEvent->GetTotalPhotons().Z();
        totCubes = dataEvent->GetTotalCubes();

        if(ind<ftrackLenghts.size() && ind<ftrackPhotos.size())
        {
            std::vector<Float_t>& tr = ftrackLenghts[ind];
            std::vector<Float_t>& ph_tr = ftrackPhotos[ind];
    
            tr1 = tr[0];
            ph_tr1 = ph_tr[0];

            tr2 = tr[1];
            ph_tr2 = ph_tr[1];

            tr3 = tr[2];
            ph_tr3 = ph_tr[2];
        }

        // if(tr1==0)
        //     continue;

        // Normalize before filling
        // tr1 /= fMaxtrack;
        // tr2 /= fMaxtrack;
        // tr3 /= fMaxtrack;
        // ph_tr1 /= fMaxTrph;
        // ph_tr2 /= fMaxTrph;
        // ph_tr3 /= fMaxTrph;
        // totCubes /= fMaxCubes;

        // totPh /= fMaxTotph;   
        // Float_t halfMaxph = fMaxTotph/2;
        // totPh = (totPh - halfMaxph)/halfMaxph; // normalize to [-1:1]


        // lnuEnergy /= fMaxnuE;
        // totalEdep /= fMaxTotEdep;

        //trueEdep /= fMaxTrueEdep
        // Float_t halfMaxTrueEdep = fMaxTrueEdep/2;
        // trueEdep = (trueEdep - halfMaxTrueEdep)/halfMaxTrueEdep; // normalize to [-1:1]

        // totalPe /= fMaxTotPe;

        longestTrackPrjTree->Fill();

        // Clear data for next write
        tr1 = 0;
        tr2 = 0;
        tr3 = 0;
        ph_tr1 = 0;
        ph_tr2 = 0;
        ph_tr3 = 0;
     }
    outFile->WriteTObject(longestTrackPrjTree);                      
    // =================================================================
    

    outFile->Close();
    
    delete outFile;

    FgdMCGenFitRecon::FinishTask();
}

Bool_t FgdTMVAData::isParticleAllowed(Int_t pdg)
{
    Bool_t isAllowed =  (pdg ==  genie::kPdgPiP) ||
                        (pdg ==  genie::kPdgPiM) ||
                      genie::pdg::IsProton(pdg) ||
                      genie::pdg::IsChargedLepton(pdg);

    return isAllowed;
}

// -------------------------------------------------------------------------


// -----   Protected methods   --------------------------------------------

// -------------------------------------------------------------------------

}// namespace superfgd
}// namespace reconstruction
}// namespace esbroot
