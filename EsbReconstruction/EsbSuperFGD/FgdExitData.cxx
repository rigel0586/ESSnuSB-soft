#include "EsbReconstruction/EsbSuperFGD/FgdExitData.h"



namespace esbroot {
namespace reconstruction {
namespace superfgd {

FgdExitData::FgdExitData()
{
}

FgdExitData::FgdExitData(Int_t nuPdg,
                Double_t nuE,
                Bool_t IsWeakCC,
                Bool_t IsWeakNC,
                Bool_t IsQuasiElastic,
                std::vector<FgdExitParticle> vecPars)
        : fnuPdg(nuPdg)
        , fnuE(nuE)
        , fIsWeakCC(IsWeakCC)
        , fIsWeakNC(IsWeakNC)
        , fIsQuasiElastic(IsQuasiElastic)
        , fVecPars(vecPars)
{
    
}

FgdExitData::~FgdExitData()
{

}

FgdExitData::FgdExitData(const FgdExitData& c)
        : fnuPdg(c.fnuPdg)
        , fnuE(c.fnuE)
        , fIsWeakCC(c.fIsWeakCC)
        , fIsWeakNC(c.fIsWeakNC)
        , fIsQuasiElastic(c.fIsQuasiElastic)
        , fVecPars(c.fVecPars)
{
    
}

FgdExitData& FgdExitData::operator=(const FgdExitData& c)
{
    fnuPdg = c.fnuPdg;
    fnuE = c.fnuE;
    fIsWeakCC = c.fIsWeakCC;
    fIsWeakNC = c.fIsWeakNC;
    fIsQuasiElastic = c.fIsQuasiElastic;
    fVecPars = c.fVecPars;
    return *this;
}

} //superfgd
} //reconstruction
} //esbroot