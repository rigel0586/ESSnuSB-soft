#ifndef ESBROOT_ESB_GRAPH_RECONSTRUCTION_FGD_GENFIT_H
#define ESBROOT_ESB_GRAPH_RECONSTRUCTION_FGD_GENFIT_H

// EsbRoot headers
#include "EsbData/EsbSuperFGD/FgdHit.h"
#include "EsbGeometry/EsbSuperFGD/EsbFgdDetectorParameters.h"
#include "EsbGeometry/EsbSuperFGD/EsbSuperFGDDetectorConstruction.h"
#include "EsbReconstruction/EsbSuperFGD/FgdReconHit.h"
#include "EsbReconstruction/EsbSuperFGD/PdgFromMomentumLoss.h"
#include "EsbReconstruction/EsbSuperFGD/FgdTMVAEventRecord.h"

// FairRoot headers
#include <FairTask.h>

// ROOT headers
#include <include/EventDisplay.h>

// Pathfinder headers
#include "basicHit.h"
#include "TrackFinderTrack.h"

namespace esbroot {
namespace reconstruction {
namespace superfgd {

class FgdGraphGenFitRecon : public FairTask
{

 public:

  enum PDG_FIT_TRACK
  {
    MUON = 13,
    PION = 211,
    PROTON = 2212,
    ELECTRON = 11
  };

  /** Default constructor **/  
  FgdGraphGenFitRecon();

  /** Constructor with argument
   *@param name       Name of task
   *@param geoConfigFile  - Configuration file detector
   *@param mediaFile  - Configuration file for the used mediums
   *@param eventData  - events data file (generated from fgd generator)
   *@param outputRootFile - full path to the output root file
   *@param verbose  - Verbosity level
   *@param debugLlv - debug level for genfit
   *@param visualize -  to visualize the event using genfit::EventDisplay
   *@param visOption -  option to be passed to genfit::EventDisplay
  **/  
  FgdGraphGenFitRecon(const char* name
              , const char* geoConfigFile
              , const char* mediaFile
              , const char* eventData
              , const char* outputRootFile
              , Int_t verbose = 1
              , double debugLlv = 0
              , bool visualize = false
              , std::string visOption ="D");

  /** Destructor **/
  virtual ~FgdGraphGenFitRecon();

  void SetMinInterations(Int_t minIterations) {fminGenFitInterations = minIterations;}
  void SetMaxInterations(Int_t maxIterations) {fmaxGenFitIterations = maxIterations;}
  void SetMinHits(Int_t minHits) {fminHits = minHits;}
  void AddPdgMomLoss(Int_t pdg, Double_t momLoss, Double_t allowDiff){ fpdgFromMomLoss.emplace_back(pdg, momLoss, allowDiff);}
  void SetSmoothCoor(bool val){fuseSmoothPos = val;}

  /** Virtual method Init **/
  virtual InitStatus Init() override;
  virtual void OutputFileInit(FairRootManager* manager); // Create virtual method for output file creation
  virtual void FinishEvent() override;
  virtual void FinishTask() override;


  /** Virtual method Exec **/
  virtual void Exec(Option_t* opt) override;

protected:

  enum FindTrackType{
    STRAIGHT_LINE,
    HELIX,
    CURL
  };

  /** Get all hits **/
  Bool_t GetHits(std::vector<ReconHit>& allHits);

  /** Extrack tracks using graph traversal and track gradient  **/
  Bool_t FindUsingGraph(std::vector<ReconHit>& hits
                  , std::vector<std::vector<ReconHit>>& foundTracks);
  
  //void BuildGraph(std::vector<ReconHit>& hits);
  void CalculateGrad(std::vector<std::vector<ReconHit*>>& tracks);
  void SplitTrack(std::vector<std::vector<ReconHit*>>& originalTracks, std::vector<std::vector<ReconHit*>>& splitTracks);
  Bool_t CalculateInitialMomentum(const std::vector<ReconHit>& track,const TVector3& magField, TVector3& momentum, TVector3& momentumLoss);
  Bool_t CalculateMomentum(const TVector3& trackVector, const TVector3& p1, const TVector3& p2, const TVector3& p3 , const TVector3& magField, TVector3& momentum);
  Bool_t CalculateCalorimetricMomentum(const std::vector<ReconHit>& track, TVector3& momentum);
  Double_t GetRadius(const TVector3& p1, const TVector3& p2, const TVector3& p3);
  Int_t GetPdgCode(const TVector3& momentum, const TVector3& momentumLoss, Int_t tryFit, Int_t& momCoeff);
  void ConvertHitToVec(std::vector<TVector3>& points, std::vector<ReconHit>& hits);

  /** Fit the found tracks using genfit **/
  void FitTracks(std::vector<std::vector<ReconHit>>& foundTracks);
  bool FitTrack(std::vector<ReconHit>& track
    , std::shared_ptr<genfit::AbsKalmanFitter>& fitter
    , Int_t pdg
    , TVector3& momentum
    , TVector3& calMom
    , int trackId);

  void ExtractTMVAdata(std::vector<ReconHit>& allhits);

  /** Define materials used in the reconstruction phase **/
  void DefineMaterials();

  /** Print information for fitted grack **/
  void PrintFitTrack(genfit::Track& track);

  Long_t ArrInd(int i, int j, int k);

  /** Class to hold the Detector parameters read from external file **/
  esbroot::geometry::superfgd::FgdDetectorParameters fParams;

  /** Class containing the TGeometry for reconstruction of the tracks **/
  esbroot::geometry::superfgd::SuperFGDDetectorConstruction    fgdConstructor;	   //! SuperFgd Detector Constructor

  /** Detector geometry pointer**/
  TGeoVolume* fsuperFgdVol;

  /** Input array of FgdDetectorPoint(s)**/
  TClonesArray* fHitArray;     //! 

  /** Output array with genfit::Track(s) **/
  TClonesArray* fTracksArray;        //!

   /** Detector dimentions **/
  Double_t flunit;
  double f_step_X;
  double f_step_Y;
  double f_step_Z;
  int f_bin_X;
  int f_bin_Y;
  int f_bin_Z;
  double f_total_X;
  double f_total_Y;
  double f_total_Z;

  /** Start position and momentum **/
  //  0 - no additional info
  //  1 - debug info
  //  2- detail info
  double fDebuglvl_genfit;

  /** Path to the used media.geo file - containing definitions of materials **/
  std::string fmediaFile;

  /** Path to the configuration file **/
  std::string fgeoConfFile;

  /** Path to the used data file containing genie events infos **/
  std::string feventData;//!<!
  std::string foutputRootFile;//!<!
  std::vector<FgdTMVAEventRecord> feventRecords;//!<!
  int feventNum;//!<!

  Int_t fminGenFitInterations;
  Int_t fmaxGenFitIterations;
  Int_t fminHits;

  /** Are materials already defined **/
  bool isDefinedMaterials;
  bool fuseSmoothPos;

  /** local Members for genfit visualization  **/
  genfit::EventDisplay* fdisplay;//!<!
  bool isGenFitVisualization;//!<!
  std::string fGenFitVisOption;//!<!

  std::vector<PdgFromMomentumLoss> fpdgFromMomLoss;//!<!

  	   
  ClassDef(FgdGraphGenFitRecon, 2);

};

} //superfgd
} //reconstruction
} //esbroot

#endif // ESBROOT_ESB_GRAPH_RECONSTRUCTION_FGD_GENFIT_H
