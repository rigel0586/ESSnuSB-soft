#ifndef ESBROOT_ESBDRECONSTRUCTION_FGD_EXITDATA_H
#define ESBROOT_ESBDRECONSTRUCTION_FGD_EXITDATA_H

#include "TObject.h"
#include "TClonesArray.h"              
#include <TVector3.h>
#include <vector>
#include <algorithm>

#include "EsbReconstruction/EsbSuperFGD/FgdExitParticle.h"

namespace esbroot {
namespace reconstruction {
namespace superfgd {

class FgdExitData : public TObject
{

public:

    /** Default constructor **/  
    FgdExitData();

    FgdExitData(Int_t nuPdg,
                Double_t nuE,
                Bool_t IsWeakCC,
                Bool_t IsWeakNC,
                Bool_t IsQuasiElastic,
                std::vector<FgdExitParticle> vecPars);

    ~FgdExitData();

    FgdExitData(const FgdExitData& c);

    FgdExitData& operator=(const FgdExitData& c);

    Int_t fnuPdg;
    Double_t fnuE;
    Bool_t fIsWeakCC;
    Bool_t fIsWeakNC;
    Bool_t fIsQuasiElastic;
    std::vector<FgdExitParticle> fVecPars;

private:

    ClassDef(FgdExitData, 2);
};

} //superfgd
} //reconstruction
} //esbroot

#endif