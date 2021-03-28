#ifndef ESBROOT_ESBDRECONSTRUCTION_FGD_MC_GENFIT_H
#define ESBROOT_ESBDRECONSTRUCTION_FGD_MC_GENFIT_H

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

class FgdMCGenFitRecon : public FairTask
{

 public:

  /** Default constructor **/  
  FgdMCGenFitRecon();

  /** Constructor with argument
   *@param name       Name of task
   *@param geoConfigFile  - Configuration file detector
   *@param mediaFile  - Configuration file for the used mediums
   *@param eventData  - events data file (generated from fgd generator)
   *@param verbose  - Verbosity level
   *@param debugLlv - debug level for genfit
   *@param visualize -  to visualize the event using genfit::EventDisplay
   *@param visOption -  option to be passed to genfit::EventDisplay
  **/  
  FgdMCGenFitRecon(const char* name
              , const char* geoConfigFile
              , const char* mediaFile
              , const char* eventData
              , Int_t verbose = 1
              , double debugLlv = 0
              , bool visualize = false
              , std::string visOption ="D");

  /** Destructor **/
  virtual ~FgdMCGenFitRecon();

  void SetMinInterations(Int_t minIterations) {fminGenFitInterations = minIterations;}
  void SetMaxInterations(Int_t maxIterations) {fmaxGenFitIterations = maxIterations;}
  void SetMinHits(Int_t minHits) {fminHits = minHits;}

  /** Virtual method Init **/
  virtual InitStatus Init() override;
  virtual void OutputFileInit(FairRootManager* manager); // Create virtual method for output file creation
  virtual void FinishEvent() override;
  virtual void FinishTask() override;


  /** Virtual method Exec **/
  virtual void Exec(Option_t* opt) override;

  /** For decendents who want to write some data to output file */
  virtual void WriteOutput( Int_t pdg
                          , const TVector3& fitMom
                          , const TVector3& mcMom
                          , const genfit::Track& fitTrack
                          , genfit::FitStatus*& fiStatuStatus
                          , std::vector<ReconHit>& track);

  void SetOutputRootFile(std::string rootFIlePath){foutputRootFile = rootFIlePath;}
  void SetErrOutFile(std::string errOutPath){ ferrOutPath = errOutPath;}
  void SetErrMomBinWidth(float momErrWidth){ fmomErrWidth = momErrWidth;}

protected:

  enum Axis
  {
    X,
    Y,
    Z
  };

  struct FitData
  {
    TVector3 fitMom;
    TVector3 mcMom;
    TVector3 trackProj;
    Int_t totalPath;
    Bool_t isFitted;

    FitData(){}
    ~FitData(){}

    FitData(const FitData& c)
    {
      fitMom = c.fitMom;
      mcMom = c.mcMom;
      trackProj = c.trackProj;
      totalPath = c.totalPath;
      isFitted = c.isFitted;
    }

    FitData& operator=(const FitData& c)
    {
      fitMom = c.fitMom;
      mcMom = c.mcMom;
      trackProj = c.trackProj;
      totalPath = c.totalPath;
      isFitted = c.isFitted;
    }
  };

  struct HoughPoint
  {
    int a;
    int b;
    int r;
    std::set<Long_t> count;

    HoughPoint(): a(0), b(0), r(0) {}
    ~HoughPoint(){}

    HoughPoint(const HoughPoint& c)
    {
      a = c.a;
      b = c.b;
      r = c.r;
      count = c.count;
    }

    HoughPoint& operator=(const HoughPoint& c)
    {
      a = c.a;
      b = c.b;
      r = c.r;
      count = c.count;
    }

    bool operator==(const HoughPoint& c)
    {
      return (a == c.a) && (b == c.b) && (r == c.r);
    }
  };

  ///      I  |  II
  ///  *------------
  ///     III | IV
  /// Starts from * and is clockwise
  /// I - from 0 to 90
  /// II - from 90 to 180
  /// III - from 180 to 270
  /// IV - from 270 to 360
  enum Quadrant
  {
      One,
      Two,
      Three,
      Four
  };

  ///         I
  ///  *------------
  ///         II
  /// Starts from * and is clockwise
  /// I - from 0 to 180
  /// II - from 180 to 360
  enum SemiCircle
  {
      SemiOne,
      SemiTwo
  };

  /** Get all hits **/
  Bool_t GetHits(std::vector<ReconHit>& allHits);

  void SplitTrack(std::vector<ReconHit>& allHits, std::vector<std::vector<ReconHit>>& splitTracks);

  /** Fit the found tracks using genfit **/
  void FitTracks(std::vector<std::vector<ReconHit>>& foundTracks);

  /** Define materials used in the reconstruction phase **/
  void DefineMaterials();

  /** Print information for fitted grack **/
  void PrintFitTrack(genfit::Track& track);

  Long_t ArrInd(int i, int j, int k);
  Long_t HitId(ReconHit& hit);
  void FillHitIds(std::vector<ReconHit>& track);

  Bool_t isParticleNeutral(Int_t pdg);
  Bool_t isAllowed(Int_t pdg);
  void   writeErrFile(const std::string& fileEnding, Axis axis, std::vector<FitData>& dataVec);
  void   writeErrFile(const std::string& fileEnding, std::vector<FitData>& dataVec);
  void   writeTrErrFile(const std::string& fileEnding, std::vector<FitData>& dataVec);
  float  calcStdDev(float& mean, std::vector<float>&);

  /** Class to hold the Detector parameters read from external file **/
  esbroot::geometry::superfgd::FgdDetectorParameters fParams;//!<!

  /** Class containing the TGeometry for reconstruction of the tracks **/
  esbroot::geometry::superfgd::SuperFGDDetectorConstruction    fgdConstructor;//!<!	   //! SuperFgd Detector Constructor

  /** Detector geometry pointer**/
  TGeoVolume* fsuperFgdVol;//!<!

  /** Input array of FgdDetectorPoint(s)**/
  TClonesArray* fHitArray;     //! 

  /** Output array with genfit::Track(s) **/
  TClonesArray* fTracksArray;        //!

   /** Detector dimentions **/
  Double_t flunit;//!<!
  double f_step_X;//!<!
  double f_step_Y;//!<!
  double f_step_Z;//!<!
  int f_bin_X;//!<!
  int f_bin_Y;//!<!
  int f_bin_Z;//!<!
  double f_total_X;//!<!
  double f_total_Y;//!<!
  double f_total_Z;//!<!

  /** Start position and momentum **/
  //  0 - no additional info
  //  1 - debug info
  //  2- detail info
  double fDebuglvl_genfit;//!<!

  /** Path to the used media.geo file - containing definitions of materials **/
  std::string fmediaFile;//!<!

  /** Path to the events file containing Monte carlo simulation data - pdg codes, momentums etc. **/
  std::string feventFile;//!<!

  Int_t fminGenFitInterations;
  Int_t fmaxGenFitIterations;
  Int_t fminHits;

  /** Are materials already defined **/
  bool isDefinedMaterials;

  /** local Members for genfit visualization  **/
  genfit::EventDisplay* fdisplay;//!<!
  bool isGenFitVisualization;//!<!
  std::string fGenFitVisOption;//!<!
  std::vector<genfit::Track*> fgenTracks;//!<!

  std::vector<Float_t> felectronFitErr;//!<!
  std::vector<Float_t> fmuonFitErr;//!<!
  std::vector<Float_t> fprotonFitErr;//!<!
  std::vector<Float_t> fpionFitErr;//!<!

  std::string fMCeventData;//!<!
  std::vector<FgdTMVAEventRecord> fMCeventRecords;//!<!
  std::string foutputRootFile;//!<!
  int fMCeventNum;//!<!

  
  std::string ferrOutPath;//!<!
  // Mom error width in MB
  float fmomErrWidth;//!<!
  std::vector<FitData> fmuonFitData;//!<!
  std::vector<FitData> fprotonFitData;//!<!
  std::vector<FitData> felectronFitData;//!<!
  std::vector<FitData> fpionFitData;//!<!
  
  	   
  ClassDef(FgdMCGenFitRecon, 2);

private:
  Double_t Radius(const TVector3& p1, const TVector3& p2, const TVector3& p3);
  bool GetHoughMomentum(std::vector<ReconHit>& track, Double_t& momentum);
  bool HoughRadius(std::vector<ReconHit>& track, Double_t& initialRadius, Double_t& radius, Int_t& cubesInCirlcle);
};

} //superfgd
} //reconstruction
} //esbroot

#endif // ESBROOT_ESBDRECONSTRUCTION_FGD_MC_GENFIT_H
