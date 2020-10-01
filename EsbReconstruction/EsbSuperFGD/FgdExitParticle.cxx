#include "EsbReconstruction/EsbSuperFGD/FgdExitParticle.h"



namespace esbroot {
namespace reconstruction {
namespace superfgd {

FgdExitParticle::FgdExitParticle()
{
}

FgdExitParticle::FgdExitParticle(Int_t Pdg,
                TVector3 Momentum,
                TVector3 position)
        : fPdg(Pdg)
        , fMomentum(Momentum)
        , fposition(position)
{
    
}

FgdExitParticle::~FgdExitParticle()
{

}

FgdExitParticle::FgdExitParticle(const FgdExitParticle& c)
        : fPdg(c.fPdg)
        , fMomentum(c.fMomentum)
        , fposition(c.fposition)
{
    
}

FgdExitParticle& FgdExitParticle::operator=(const FgdExitParticle& c)
{
    fPdg = c.fPdg;
    fMomentum = c.fMomentum;
    fposition = c.fposition;
    return *this;
}

} //superfgd
} //reconstruction
} //esbroot