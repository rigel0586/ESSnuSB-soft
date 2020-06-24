#include "EsbReconstruction/EsbSuperFGD/FgdTMVAData2.h"
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
#define MAX_PHOTONS_NORMALIZE 1000 // Number of photons to normalize

// -----   Default constructor   -------------------------------------------
FgdTMVAData2::FgdTMVAData2() : FgdMCGenFitRecon(), feventNum(0), fmagField_X(0.), fmagField_Y(0.), fmagField_Z(0.)
{
}
// -------------------------------------------------------------------------

// -----   Constructor   -------------------------------------------
FgdTMVAData2::FgdTMVAData2(const char* name
                          , const char* geoConfigFile
                          , const char* mediaFile
                          , const char* eventData
                          , const char* outputRootFile
                          , Int_t verbose
                          , double debugLlv) :
  FgdMCGenFitRecon(name, geoConfigFile, mediaFile, eventData, verbose, 
                    debugLlv, false /* no visualization */, "D")
    , feventData(eventData), foutputRootFile(outputRootFile)
    , feventNum(0), fmagField_X(0.), fmagField_Y(0.), fmagField_Z(0.)
    ,  foutFile(nullptr), fdeTree(nullptr), fciTree(nullptr)
{ 
}
// -------------------------------------------------------------------------



// -----   Destructor   ----------------------------------------------------
FgdTMVAData2::~FgdTMVAData2() 
{
    if(fTracksArray) 
    {
        fTracksArray->Delete();
        delete fTracksArray;
    } 

    if(foutFile!=nullptr)
    {
        delete foutFile;
    }
}
// -------------------------------------------------------------------------



// -----   Public method Init   --------------------------------------------
InitStatus FgdTMVAData2::Init() 
{   
    FgdMCGenFitRecon::Init();

    fmagField_X = fParams.ParamAsDouble(esbroot::geometry::superfgd::DP::magField_X);
    fmagField_Y = fParams.ParamAsDouble(esbroot::geometry::superfgd::DP::magField_Y);
    fmagField_Z = fParams.ParamAsDouble(esbroot::geometry::superfgd::DP::magField_Z);

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

    foutFile = new TFile(foutputRootFile.c_str(), "RECREATE", "TVMA data from Fgd Detector");

    fdeTree = new TTree(esbroot::geometry::superfgd::DP::FGD_DEPOSITED_ENERGY_CUBE_TTREE.c_str()
                                ,esbroot::geometry::superfgd::DP::FGD_TMVA_DATA_ROOT_FILE.c_str());

    fciTree = new TTree(esbroot::geometry::superfgd::DP::FGD_CUBE_INFO_TTREE.c_str()
                                ,esbroot::geometry::superfgd::DP::FGD_TMVA_DATA_ROOT_FILE.c_str());

    return kSUCCESS;
}

void FgdTMVAData2::OutputFileInit(FairRootManager* manager)
{
    FgdMCGenFitRecon::OutputFileInit(manager);
}

// -------------------------------------------------------------------------



// -----   Public methods   --------------------------------------------

void FgdTMVAData2::Exec(Option_t* opt) 
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

Bool_t FgdTMVAData2::ProcessStats(std::vector<std::vector<ReconHit>>& foundTracks)
{
    if(feventNum >= feventRecords.size())
    {
        LOG(fatal) << "EventData reconrds are less than the passed events!";
        throw "EventData reconrds are less than the passed events!";
    }

    FgdTMVAEventRecord& tvmaEventRecord = feventRecords[feventNum];
    tvmaEventRecord.ReadEventData();

    if(!tvmaEventRecord.IsWeakCC() || !tvmaEventRecord.IsQuasiElastic())
    {
        ++feventNum;
        return true;
    }

    Float_t dE(0);
    Float_t phot(0);
    fdeTree->Branch("dE", &dE);
    fdeTree->Branch("phot", &phot);

    Float_t photons(0);
    Float_t dedx(0);
    fdeTree->Branch("photons", &photons);
    fdeTree->Branch("dedx", &dedx);

    Float_t dst_X(0);
    Float_t dst_X2(0);
    fdeTree->Branch("dst_X", &dst_X);
    fdeTree->Branch("dst_X2", &dst_X2);

    Float_t dst_Y(0);
    Float_t dst_Y2(0);
    fdeTree->Branch("dst_Y", &dst_Y);
    fdeTree->Branch("dst_Y2", &dst_Y2);

    Float_t dst_Z(0);
    Float_t dst_Z2(0);
    fdeTree->Branch("dst_Z", &dst_Z);
    fdeTree->Branch("dst_Z2", &dst_Z2);

    Int_t totalCubes(0);
    fciTree->Branch("totalCubes", &totalCubes);
    

    for(size_t i = 0; i <  foundTracks.size() ; ++i)
    {
        std::vector<ReconHit>& hitsOnTrack = foundTracks[i];
        if(hitsOnTrack.empty()) continue;

        for(size_t j = 0; j < hitsOnTrack.size(); ++j)
        {   
            ReconHit& hit = hitsOnTrack[j];

            photons = hit.fpe;
            dedx = hit.fEdep / hit.ftrackLength;

            dE = hit.fEdep;
            phot =   hit.fph1.X() + hit.fph1.Y() + hit.fph1.Z()
                   + hit.fph2.X() + hit.fph2.Y() + hit.fph2.Z();

            phot /= MAX_PHOTONS_NORMALIZE;

            TVector3& dist_1 = hit.fmppc1;
            dst_X = dist_1.X()/f_bin_X;
            dst_Y = dist_1.Y()/f_bin_Y;
            dst_Z = dist_1.Z()/f_bin_Z;

            TVector3& dist_2 = hit.fmppc2;
            dst_X2 = (1.*dist_2.X())/f_bin_X;
            dst_Y2 = (1*dist_2.Y())/f_bin_Y;
            dst_Z2 = (1*dist_2.Z())/f_bin_Z;

            fdeTree->Fill();
            ++totalCubes;
        }   
    }

    fciTree->Fill();
    ++feventNum; // Increment to next event from eventData read from simulation`s genie export
}

void FgdTMVAData2::FinishTask()
{
    foutFile->WriteTObject(fdeTree); 
    foutFile->WriteTObject(fciTree); 
    FgdMCGenFitRecon::FinishTask();
}

// -------------------------------------------------------------------------


// -----   Protected methods   --------------------------------------------

// -------------------------------------------------------------------------

}// namespace superfgd
}// namespace reconstruction
}// namespace esbroot
