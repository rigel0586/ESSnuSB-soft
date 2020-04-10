#include "EsbReconstruction/EsbSuperFGD/FgdEdepAnalyzer.h"
#include "EsbReconstruction/EsbSuperFGD/FgdReconTemplate.h"
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
FgdEdepAnalyzer::FgdEdepAnalyzer() : FgdMCGenFitRecon(), feventNum(0)
        , fmagField_X(0.), fmagField_Y(0.), fmagField_Z(0.), fDataArr()
{
}
// -------------------------------------------------------------------------

// -----   Constructor   -------------------------------------------
FgdEdepAnalyzer::FgdEdepAnalyzer(const char* name
                          , const char* geoConfigFile
                          , const char* mediaFile
                          , const char* outputEdepFile
                          , Int_t photoInterval
                          , Int_t verbose
                          , double debugLlv) :
  FgdMCGenFitRecon(name, geoConfigFile, mediaFile, verbose, 
                    debugLlv, false /* no visualization */, "D")
    , foutputEdepFile(outputEdepFile), fphotoInterval(photoInterval)
    , feventNum(0), fmagField_X(0.), fmagField_Y(0.), fmagField_Z(0.)
    , fDataArr()
{ 
    fphotoInterval = (fphotoInterval <=0)? 1: fphotoInterval;
    fDataArr.setPhInt(fphotoInterval);
}
// -------------------------------------------------------------------------



// -----   Destructor   ----------------------------------------------------
FgdEdepAnalyzer::~FgdEdepAnalyzer() 
{
    if(fTracksArray) 
    {
        fTracksArray->Delete();
        delete fTracksArray;
    } 
}
// -------------------------------------------------------------------------



// -----   Public method Init   --------------------------------------------
InitStatus FgdEdepAnalyzer::Init() 
{   
    FgdMCGenFitRecon::Init();

    fmagField_X = fParams.ParamAsDouble(esbroot::geometry::superfgd::DP::magField_X);
    fmagField_Y = fParams.ParamAsDouble(esbroot::geometry::superfgd::DP::magField_Y);
    fmagField_Z = fParams.ParamAsDouble(esbroot::geometry::superfgd::DP::magField_Z); 

    return kSUCCESS;
}

void FgdEdepAnalyzer::OutputFileInit(FairRootManager* manager)
{
    FgdMCGenFitRecon::OutputFileInit(manager);
}

// -------------------------------------------------------------------------



// -----   Public methods   --------------------------------------------

void FgdEdepAnalyzer::Exec(Option_t* opt) 
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
    }
  }
  catch(genfit::Exception& e)
  {
      LOG(fatal) << "Exception, when tryng to fit track";
      LOG(fatal) << e.what();
  }
}

Bool_t FgdEdepAnalyzer::ProcessStats(std::vector<std::vector<ReconHit>>& foundTracks)
{
    for(size_t i = 0; i <  foundTracks.size() ; ++i)
    {
        std::vector<ReconHit>& hitsOnTrack = foundTracks[i];
        if(hitsOnTrack.empty()) continue;

        for(size_t j = 0; j < hitsOnTrack.size(); ++j)
        {   
            ReconHit& hit = hitsOnTrack[j];
            Double_t trackLength = hit.ftrackLength;
            Double_t edep = hit.fEdep;
            Double_t pe = hit.fpe;

            //LOG(warning) << "pe " << pe << " dedx " <<  edep/trackLength << " [ " << edep << ", " << trackLength << "]"; 
            fDataArr.add(pe, edep/trackLength);
        }   
    }
}

void FgdEdepAnalyzer::FinishTask()
{
    FgdMCGenFitRecon::FinishTask();

    LOG(warning) << "==========================================";
    std::vector<EdepInfo>& info = fDataArr.getInfo();
    std::sort(info.begin(), info.end(), [](EdepInfo& bh1, EdepInfo& bh2){return bh1.getLow()<bh2.getLow();});
    for(size_t i = 0; i <  info.size() ; ++i)
    {
        EdepInfo& elem = info[i];
        LOG(warning) << "[ " <<  elem.getLow() << "," << elem.getHigh() << "," << elem.size() << "] = " << elem.getAvgDedx();
    }
}

FgdEdepAnalyzer::EdepInfo::EdepInfo(Int_t low, Int_t high)
    : flow(low), fhigh(high)
{

}

FgdEdepAnalyzer::EdepInfo::EdepInfo(const EdepInfo& c)
{
    flow = c.flow;
    fhigh = c.fhigh;
    fdedx = c.fdedx;
}

FgdEdepAnalyzer::EdepInfo::~EdepInfo()
{

}

Bool_t FgdEdepAnalyzer::EdepInfo::contains(Double_t pe)
{
    return (flow<= pe && pe < fhigh);
}

Double_t FgdEdepAnalyzer::EdepInfo::getAvgDedx()
{
    Double_t sum = 0.;
    for(size_t i =0; i < fdedx.size(); ++i)
    {
        sum+= fdedx[i];
    }
    Double_t&& avg = sum / fdedx.size();
    return avg;
}


void FgdEdepAnalyzer::EdepArray::add(Double_t pe, Double_t dedx)
{
    Bool_t containsVal(false);
    for(size_t i =0; i < fInfo.size(); ++i)
    {
        if(fInfo[i].contains(pe))
        {
            fInfo[i].add(dedx);
            containsVal = true;
            break;
        }
    }

    if(!containsVal)
    {
        Int_t low = ((Int_t)pe/fphInt) * fphInt;
        Int_t high = low + fphInt;
        EdepInfo elem(low , high);
        elem.add(dedx);
        fInfo.emplace_back(elem);
    }
}

// -------------------------------------------------------------------------


// -----   Protected methods   --------------------------------------------

// -------------------------------------------------------------------------

}// namespace superfgd
}// namespace reconstruction
}// namespace esbroot
