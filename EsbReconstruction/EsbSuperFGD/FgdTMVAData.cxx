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

// -----   Default constructor   -------------------------------------------
FgdTMVAData::FgdTMVAData() : FgdMCGenFitRecon(), feventNum(0)
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
    , feventNum(0)
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

    for(size_t i = 0; i <  foundTracks.size() ; ++i)
    {
        std::vector<ReconHit>& hitsOnTrack = foundTracks[i];
        // Sort by time, the 1st hit in time is the start of the track
        std::sort(hitsOnTrack.begin(), hitsOnTrack.end(), [](ReconHit& bh1, ReconHit& bh2){return bh1.ftime<bh2.ftime;});
    }
    
    for(size_t i = 0; i <  foundTracks.size() ; ++i)
    {
        std::vector<ReconHit>& hitsOnTrack = foundTracks[i];
        if(hitsOnTrack.empty()) continue;

        Int_t sumTotalPhoto = 0;
        Int_t sumTotalCubes = 0;
        for(size_t j = 0; j < hitsOnTrack.size(); ++j)
        {   
            ReconHit& hit = hitsOnTrack[j];
            sumTotalPhoto += hit.fphotons.X() + hit.fphotons.Y() + hit.fphotons.Z();
            sumTotalCubes++;

            fhitCoordinates[feventNum].emplace_back(hit.fHitPos);
            fhitPhotons[feventNum].emplace_back(hit.fphotons);
        }

        tvmaEventRecord.SetTotalPhotons(sumTotalPhoto);
        tvmaEventRecord.SetTotalCubes(sumTotalCubes);
    }

    feventRecords[feventNum].SetHasHits(true);
    ++feventNum; // Increment to next event from eventData read from simulation`s genie export
}



void FgdTMVAData::FinishTask()
{
    TFile * outFile = new TFile(foutputRootFile.c_str(), "RECREATE", "TVMA data from Fgd Detector");
    outFile->SetCompressionLevel(9);

    FgdTMVAEventRecord* data = nullptr;
    TClonesArray* hitCoordinates = new TClonesArray(TVector3::Class());
    TClonesArray& hitcref = *hitCoordinates;

    TClonesArray* hitPhotons = new TClonesArray(TVector3::Class());
    TClonesArray& photoref = *hitPhotons;

    TTree * outTree = new TTree(esbroot::geometry::superfgd::DP::FGD_TMVA_DATA_TTREE.c_str()
                                ,esbroot::geometry::superfgd::DP::FGD_TMVA_DATA_ROOT_FILE.c_str());
    outTree->Branch(esbroot::geometry::superfgd::DP::FGD_TMVA_DATA_BRANCH.c_str(), &data);
    outTree->Branch(esbroot::geometry::superfgd::DP::FGD_TMVA_HIT_ARRAY_BRANCH.c_str(), &hitCoordinates);
    outTree->Branch(esbroot::geometry::superfgd::DP::FGD_TMVA_PHOTO_ARRAY_BRANCH.c_str(), &hitPhotons);

 
    for(size_t ind = 0 ; ind < feventRecords.size(); ind++)
    {
        data = &feventRecords[ind];

        // 1. Copy hitposition
        for(Int_t i = 0 ; i < fhitCoordinates[ind].size(); i++)
        {
            new(hitcref[i]) TVector3(fhitCoordinates[ind][i]);
        }

        // 2. Copy photons
        for(Int_t i = 0 ; i < fhitPhotons[ind].size(); i++)
        {
            new(photoref[i]) TVector3(fhitPhotons[ind][i]);
        }

        outTree->Fill();

        hitCoordinates->Clear();
    }

    outFile->WriteTObject(outTree);
    outFile->Close();
    
    delete outFile;

    FgdMCGenFitRecon::FinishTask();
}
// -------------------------------------------------------------------------


// -----   Protected methods   --------------------------------------------

// -------------------------------------------------------------------------

}// namespace superfgd
}// namespace reconstruction
}// namespace esbroot
