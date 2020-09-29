#ifndef ESBROOT_ESBDRECONSTRUCTION_FGD_MC_LEPTON_STATS_H
#define ESBROOT_ESBDRECONSTRUCTION_FGD_MC_LEPTON_STATS_H

// EsbRoot headers
#include "EsbReconstruction/EsbSuperFGD/FgdMCGenFitRecon.h"
#include "EsbReconstruction/EsbSuperFGD/FgdMCEventRecord.h"
#include <fstream>

namespace esbroot {
namespace reconstruction {
namespace superfgd {

class FgdMCLeptonStats : public FgdMCGenFitRecon
{

 public:

  /** Default constructor **/  
  FgdMCLeptonStats();

  /** Constructor with argument
   *@param name       Name of task
   *@param geoConfigFile  - Configuration file detector
   *@param mediaFile  - Configuration file for the used mediums
   *@param eventData  - events data file (generated from fgd generator)
   *@param outputRootFile - full path to the output root file
   *@param exitParticlesFile - output file with data for exiting particles
   *@param verbose  - Verbosity level
   *@param debugLlv - debug level for genfit
  **/  
  FgdMCLeptonStats(const char* name
              , const char* geoConfigFile
              , const char* mediaFile
              , const char* eventData
              , const char* outputRootFile
              , const char* exitParticlesFile
              , Int_t verbose = 1
              , double debugLlv = 0);

  /** Destructor **/
  virtual ~FgdMCLeptonStats();

  /** Virtual method Init **/
  virtual InitStatus Init() override;
  virtual void OutputFileInit(FairRootManager* manager) override;

  /** Virtual method Exec **/
  virtual void Exec(Option_t* opt) override;

  virtual void FinishTask() override;

  void SetFitMuonMom(bool val){ffitMuonMom=val;}

protected:

  Bool_t ProcessStats(std::vector<std::vector<ReconHit>>& foundTracks);

  Bool_t IsHitExiting(ReconHit& hit);
  Bool_t IsChargedHadron(ReconHit& hit);
  Bool_t IsChargedParticle(ReconHit& hit);
  Bool_t FitTrack(std::vector<ReconHit>& hitsOnTrack, TVector3& fitMom);
  Double_t CalculatePhotoEdep(ReconHit& hit);
  
  std::string feventData;//!<!
  std::string foutputRootFile;//!<!
  std::string fexitParticlesFile;//!<!
  std::ofstream foutputFile;//!<!
  Bool_t ffitMuonMom;//!<!

  std::vector<FgdMCEventRecord> feventRecords;//!<!

  TClonesArray* fEventsArray;     //! 

  std::shared_ptr<TDatabasePDG> fpdgDB;//!<!

  int feventNum;//!<!
  	   
  ClassDef(FgdMCLeptonStats, 2);

};

} //superfgd
} //reconstruction
} //esbroot

#endif // ESBROOT_ESBDRECONSTRUCTION_FGD_MC_LEPTON_STATS_H
