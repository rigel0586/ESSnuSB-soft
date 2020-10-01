#ifndef ESBROOT_ESBDRECONSTRUCTION_FGD_EXITPARTICLE_H
#define ESBROOT_ESBDRECONSTRUCTION_FGD_EXITPARTICLE_H

#include "TObject.h"
#include "TClonesArray.h"              
#include <TVector3.h>
#include <vector>
#include <algorithm>

namespace esbroot {
namespace reconstruction {
namespace superfgd {

class FgdExitParticle : public TObject
{

public:

    /** Default constructor **/  
    FgdExitParticle();

    FgdExitParticle(Int_t fPdg,
                TVector3 momentum,
                TVector3 position);

    ~FgdExitParticle();

    FgdExitParticle(const FgdExitParticle& c);

    FgdExitParticle& operator=(const FgdExitParticle& c);

    Int_t fPdg;
    TVector3 fMomentum;
    TVector3 fposition;

private:

    ClassDef(FgdExitParticle, 2);
};

} //superfgd
} //reconstruction
} //esbroot

#endif