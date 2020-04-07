#ifndef ESBROOT_ESBDRECONSTRUCTION_FGD_TMVA_DATA2_STATS_H
#define ESBROOT_ESBDRECONSTRUCTION_FGD_TMVA_DATA2_STATS_H

// EsbRoot headers
#include "EsbReconstruction/EsbSuperFGD/FgdMCGenFitRecon.h"
#include "EsbReconstruction/EsbSuperFGD/FgdTMVAEventRecord.h"

namespace esbroot {
namespace reconstruction {
namespace superfgd {

class FgdTMVAData2 : public FgdMCGenFitRecon
{

 public:

  /** Default constructor **/  
  FgdTMVAData2();

  /** Constructor with argument
   *@param name       Name of task
   *@param geoConfigFile  - Configuration file detector
   *@param mediaFile  - Configuration file for the used mediums
   *@param eventData  - events data file (generated from fgd generator)
   *@param outputRootFile - full path to the output root file
   *@param verbose  - Verbosity level
   *@param debugLlv - debug level for genfit
  **/  
  FgdTMVAData2(const char* name
              , const char* geoConfigFile
              , const char* mediaFile
              , const char* eventData
              , const char* outputRootFile
              , Int_t verbose = 1
              , double debugLlv = 0);

  /** Destructor **/
  virtual ~FgdTMVAData2();

  /** Virtual method Init **/
  virtual InitStatus Init() override;
  virtual void OutputFileInit(FairRootManager* manager) override;

  /** Virtual method Exec **/
  virtual void Exec(Option_t* opt) override;

  virtual void FinishTask() override;

protected:

  Bool_t ProcessStats(std::vector<std::vector<ReconHit>>& foundTracks);
  
  std::string feventData;//!<!
  std::string foutputRootFile;//!<!

  std::vector<FgdTMVAEventRecord> feventRecords;//!<!

  int feventNum;//!<!
  Float_t fmagField_X;
  Float_t fmagField_Y;
  Float_t fmagField_Z;

  TFile * foutFile;
  TTree * fdeTree;
  TTree * fciTree;
  	   
  ClassDef(FgdTMVAData2, 2);

};

} //superfgd
} //reconstruction
} //esbroot

#endif // ESBROOT_ESBDRECONSTRUCTION_FGD_TMVA_DATA2_STATS_H
