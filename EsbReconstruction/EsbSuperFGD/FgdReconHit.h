#ifndef ESBROOT_ESBDRECONSTRUCTION_FGD_RECONHIT_H
#define ESBROOT_ESBDRECONSTRUCTION_FGD_RECONHIT_H

#include "TObject.h"
#include <TVector3.h>
#include <vector>
#include <algorithm>

namespace esbroot {
namespace reconstruction {
namespace superfgd {

class ReconHit : public TObject
{

public:

    /** Default constructor **/  
    ReconHit();

    ReconHit(TVector3 mppcLoc
            , TVector3 hitPosition
            , TVector3 mCPos
            , TVector3 photons
            , Double_t time
            , TVector3 mom
            , TVector3 momExit
            , Double_t trackLength
            , Double_t trackLengthOrigin
            , Int_t pdg
            , Int_t trackId
            , Double_t edep
            , TVector3 ph1
            , TVector3 mppc1
            , TVector3 ph2
            , TVector3 mppc2
            , Double_t pe);

    ~ReconHit();

    ReconHit(const ReconHit& c);

    ReconHit& operator=(const ReconHit& c);

    Bool_t operator==(const ReconHit& c);

    TVector3 operator-(const ReconHit& c);

    Bool_t IsAlone()
    {
        return fAllHits.empty(); 
    } 


    TVector3 fmppcLoc; // from 0 to N (mppcLoc)
    TVector3 fHitPos;  // -f_total_X/2 + f_step_X*mppcLoc.X()  +f_step_X/2;
    TVector3 fMCPos;  // MC real position of the hit;
    TVector3 fphotons;
    TVector3 fph1;
    TVector3 fmppc1;
    TVector3 fph2;
    TVector3 fmppc2;
    Double_t ftime;
    TVector3 fmom;
    TVector3 fmomExit;
    Double_t ftrackLength;       // Lenght used bu summing on each step (it is less than ftrackLengthOrigin)
    Double_t ftrackLengthOrigin; // TVirtualMC::GetMC()->TrackLength() -> Monte Carlo length
    Double_t fEdep;
    Int_t fpdg;
    Int_t ftrackId;
    Bool_t fIsLeaf;

    Double_t fXaxisAngle;
    Double_t fYaxisAngle;
    Double_t fZaxisAngle;
    Double_t fChangeAngle;
    Double_t fpe;

    Int_t fLocalId;
    Bool_t fIsVisited;
    std::vector<ReconHit*> fAllHits;//!<! 

private:

    ClassDef(ReconHit, 2);
};

} //superfgd
} //reconstruction
} //esbroot

#endif