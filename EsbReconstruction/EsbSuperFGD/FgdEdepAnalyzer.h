#ifndef ESBROOT_ESBDRECONSTRUCTION_FGD_EDEP_ANALYZER_H
#define ESBROOT_ESBDRECONSTRUCTION_FGD_EDEP_ANALYZER_H

// EsbRoot headers
#include "EsbReconstruction/EsbSuperFGD/FgdMCGenFitRecon.h"
#include "EsbReconstruction/EsbSuperFGD/FgdTMVAEventRecord.h"

namespace esbroot {
namespace reconstruction {
namespace superfgd {

class FgdEdepAnalyzer : public FgdMCGenFitRecon
{

 public:

  /** Default constructor **/  
  FgdEdepAnalyzer();

  /** Constructor with argument
   *@param name       Name of task
   *@param geoConfigFile  - Configuration file detector
   *@param mediaFile  - Configuration file for the used mediums
   *@param outputEdepFile  - output data statistics for photons to edep
   *@param photoInterval - interval for which to make statistics
   *@param verbose  - Verbosity level
   *@param debugLlv - debug level for genfit
  **/  
  FgdEdepAnalyzer(const char* name
              , const char* geoConfigFile
              , const char* mediaFile
              , const char* outputEdepFile
              , Int_t photoInterval
              , Int_t verbose = 1
              , double debugLlv = 0);

  /** Destructor **/
  virtual ~FgdEdepAnalyzer();

  /** Virtual method Init **/
  virtual InitStatus Init() override;
  virtual void OutputFileInit(FairRootManager* manager) override;

  /** Virtual method Exec **/
  virtual void Exec(Option_t* opt) override;

  virtual void FinishTask() override;

protected:

  class EdepInfo{
    public:
        EdepInfo(Int_t low, Int_t high);
        EdepInfo(const EdepInfo& c);
        ~EdepInfo();


    Bool_t contains(Double_t pe);
    void add(Double_t dedx) {fdedx.push_back(dedx);};
    Int_t getLow(){return flow;}
    Int_t getHigh(){return fhigh;}
    size_t size(){return fdedx.size();}
    Double_t getAvgDedx();

    private:
        Int_t flow;
        Int_t fhigh;
        std::vector<Double_t> fdedx;
  };

  class EdepArray
  {
    public:
        EdepArray(){};
        EdepArray(Int_t phInt) : fphInt(phInt){};
        ~EdepArray(){};

        void add(Double_t pe, Double_t dedx);
        std::vector<EdepInfo>& getInfo() {return fInfo;}
        void setPhInt(Int_t val) {fphInt = val;}
    private:
        std::vector<EdepInfo> fInfo;
        Int_t fphInt;
  };

  Bool_t ProcessStats(std::vector<std::vector<ReconHit>>& foundTracks);
  
  std::string foutputEdepFile;//!<!

  int feventNum;//!<!
  Float_t fmagField_X;
  Float_t fmagField_Y;
  Float_t fmagField_Z;

  Int_t fphotoInterval;
  EdepArray fDataArr;
  	   
  ClassDef(FgdEdepAnalyzer, 2);

};

} //superfgd
} //reconstruction
} //esbroot

#endif // ESBROOT_ESBDRECONSTRUCTION_FGD_EDEP_ANALYZER_H
