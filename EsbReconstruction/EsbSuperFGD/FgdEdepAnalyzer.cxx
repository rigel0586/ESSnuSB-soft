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
        , fmagField_X(0.), fmagField_Y(0.), fmagField_Z(0.), fmapEdep()
{
}
// -------------------------------------------------------------------------

// -----   Constructor   -------------------------------------------
FgdEdepAnalyzer::FgdEdepAnalyzer(const char* name
                          , const char* geoConfigFile
                          , const char* mediaFile
                          , const char* eventData
                          , const char* outputEdepFile
                          , Int_t photoInterval
                          , Int_t verbose
                          , double debugLlv) :
  FgdMCGenFitRecon(name, geoConfigFile, mediaFile, verbose, 
                    debugLlv, false /* no visualization */, "D")
    , feventData(eventData)
    , foutputEdepFile(outputEdepFile), fphotoInterval(photoInterval)
    , feventNum(0), fmagField_X(0.), fmagField_Y(0.), fmagField_Z(0.)
    , fmapEdep()
{ 
    fphotoInterval = (fphotoInterval <=0)? 1: fphotoInterval;
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
    std::ifstream eventFileStream;
    try
    {        
        eventFileStream.open(feventData.c_str(), std::ios::in);

        if(eventFileStream.is_open())
        {
            std::string line;
            while(std::getline(eventFileStream,line))
            {
                feventRecords.emplace_back(FgdMCEventRecord(line));
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
      ++feventNum;
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
    if(feventNum >= feventRecords.size())
    {
        LOG(fatal) << "EventData reconrds are less than the passed events!";
        throw "EventData reconrds are less than the passed events!";
    }

    FgdMCEventRecord& mcEventRecord = feventRecords[feventNum];
    mcEventRecord.ReadEventData();

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

        for(size_t j = 0; j < hitsOnTrack.size(); ++j)
        {   
            ReconHit& hit = hitsOnTrack[j];
            Double_t trackLength = hit.ftrackLength;
            Double_t edep = hit.fEdep *1e3;
            Double_t pe = hit.fpe;

            //LOG(warning) << "pdg " << hit.fpdg << " pe " << pe << " dedx " <<  edep/trackLength << " [ " << edep << ", " << trackLength << "]"; 

            // LOG(warning) << "pdg " << hit.fpdg << " pe " << pe << " Digitized " 
            //             << ApplyScintiResponse(edep, trackLength, 1);  
            
            // LOG(warning)<< "pdg " << hit.fpdg << " pe " << pe << " dedx " <<  edep/trackLength << " " << hit.fmppcLoc.X()  
            //     << " " << hit.fmppcLoc.Y() 
            //     << " " << hit.fmppcLoc.Z() ;

            std::map<Int_t, EdepArray>::iterator pdgKey= fmapEdep.find(hit.fpdg);
            if(pdgKey != fmapEdep.end())
            {
                pdgKey->second.add(pe, edep/trackLength, edep, trackLength);
            }else{
                fmapEdep[hit.fpdg].setPhInt(fphotoInterval);
                fmapEdep[hit.fpdg].add(pe, edep/trackLength, edep, trackLength);
            }
        }   
    }

    feventRecords[feventNum].SetHasHits(true);
    ++feventNum; // Increment to next event from eventData read from simulation`s genie export
}

void FgdEdepAnalyzer::FinishTask()
{
    FgdMCGenFitRecon::FinishTask();

    LOG(warning) << "==========================================";
    std::map<Int_t, EdepArray>::iterator pdgKey= fmapEdep.begin();

    while(pdgKey != fmapEdep.end())
    {
        std::vector<EdepInfo>& info = fmapEdep[pdgKey->first].getInfo();
        std::sort(info.begin(), info.end(), [](EdepInfo& bh1, EdepInfo& bh2){return bh1.getLow()<bh2.getLow();});
        for(size_t i = 0; i <  info.size() ; ++i)
        {
            EdepInfo& elem = info[i];
            LOG(warning)    << pdgKey->first << " [ " <<  elem.getLow() << "," << elem.getHigh() << "," << elem.size() << "] = " 
                            << elem.getAvgDedx()
                            << " [ StDev = " << elem.getStdDev() << " %" << 100*elem.getStdDev()/ elem.getAvgDedx() << " ] "
                            << " Edep -> "
                            << elem.getAvgEdep()
                            << " [ StDev = " << elem.getStdDevEdep() << " %" << 100*elem.getStdDevEdep()/ elem.getAvgEdep() <<  " ] "
                            << " Length -> "
                            << elem.getAvgLength()
                            << " [ StDev = " << elem.getStdDevLength() << " %" << 100*elem.getStdDevLength()/ elem.getAvgLength() <<  " ] ";
        }
        ++pdgKey;
    }
    
}

double FgdEdepAnalyzer::ApplyScintiResponse(double edep, double trackLength, double charge)
{
    const double CBIRKS = 0.0208; // used in SciBooNE MC
    const double EdepToPhotConv_FGD = 70.8 / CLHEP::MeV; // contains both collection in fiber and edep->gamma conversion 
    const double EdepToPhotConv_SuperFGD = EdepToPhotConv_FGD * 1.3;

    // Calculated as in 
    //  geant4/examples/advanced/amsEcal/src/SteppingAction.cc
    //  G4double SteppingAction::BirksAttenuation(const G4Step* aStep)
    //  light attenuation (Birks' formula)
    double dedx = edep/trackLength;
    double edep_birk = edep/(1. + CBIRKS*dedx);
    // Calculate edep to photons
    LOG(info) << "dedx " << dedx 
                << " edep " << edep 
                << " trackLength " << trackLength 
                << " (1. + CBIRKS*dedx) " << (1. + CBIRKS*dedx)
                << " edep_birk " << edep_birk
                << " CLHEP::cm " << CLHEP::cm
                << " CLHEP::MeV " << CLHEP::MeV
                << " CLHEP::cm/CLHEP::MeV " << CLHEP::cm/CLHEP::MeV;
    return edep_birk * EdepToPhotConv_SuperFGD;  
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
    fEdep = c.fEdep;
    fLength = c.fLength;
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

Double_t FgdEdepAnalyzer::EdepInfo::getStdDev()
{
    Double_t mean = this->getAvgDedx();

    Double_t sq = 0.;
    for(size_t i =0; i < fdedx.size(); ++i)
    {
        sq+= std::pow( (fdedx[i] - mean), 2);
    }
    Double_t&& avgSq = sq / fdedx.size();
    return std::sqrt(avgSq);
}

Double_t FgdEdepAnalyzer::EdepInfo::getAvgEdep()
{
    Double_t sum = 0.;
    for(size_t i =0; i < fEdep.size(); ++i)
    {
        sum+= fEdep[i];
    }
    Double_t&& avg = sum / fEdep.size();
    return avg;
}

Double_t FgdEdepAnalyzer::EdepInfo::getStdDevEdep()
{
    Double_t mean = this->getAvgEdep();

    Double_t sq = 0.;
    for(size_t i =0; i < fEdep.size(); ++i)
    {
        sq+= std::pow( (fEdep[i] - mean), 2);
    }
    Double_t&& avgSq = sq / fEdep.size();
    return std::sqrt(avgSq);
}

Double_t FgdEdepAnalyzer::EdepInfo::getAvgLength()
{
    Double_t sum = 0.;
    for(size_t i =0; i < fLength.size(); ++i)
    {
        sum+= fLength[i];
    }
    Double_t&& avg = sum / fLength.size();
    return avg;
}

Double_t FgdEdepAnalyzer::EdepInfo::getStdDevLength()
{
    Double_t mean = this->getAvgLength();

    Double_t sq = 0.;
    for(size_t i =0; i < fLength.size(); ++i)
    {
        sq+= std::pow( (fLength[i] - mean), 2);
    }
    Double_t&& avgSq = sq / fLength.size();
    return std::sqrt(avgSq);
}


void FgdEdepAnalyzer::EdepInfo::add(Double_t dedx, Double_t edep, Double_t length)
{
    fdedx.push_back(dedx);
    fEdep.push_back(edep);
    fLength.push_back(length);
}


void FgdEdepAnalyzer::EdepArray::add(Double_t pe, Double_t dedx, Double_t edep, Double_t length)
{
    Bool_t containsVal(false);
    for(size_t i =0; i < fInfo.size(); ++i)
    {
        if(fInfo[i].contains(pe))
        {
            fInfo[i].add(dedx, edep, length);
            containsVal = true;
            break;
        }
    }

    if(!containsVal)
    {
        Int_t low = ((Int_t)pe/fphInt) * fphInt;
        Int_t high = low + fphInt;
        EdepInfo elem(low , high);
        elem.add(dedx, edep, length);
        fInfo.emplace_back(elem);
    }
}

// -------------------------------------------------------------------------


// -----   Protected methods   --------------------------------------------

// -------------------------------------------------------------------------

}// namespace superfgd
}// namespace reconstruction
}// namespace esbroot
