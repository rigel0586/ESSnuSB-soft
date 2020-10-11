#include "EsbReconstruction/EsbSuperFGD/FgdGraphRecord.h"



namespace esbroot {
namespace reconstruction {
namespace superfgd {

FgdGraphRecord::FgdGraphRecord()
    :     fvecDist(0)
        , fdisDiff(0)
        , fgradDiffAlowed(0.)
{
}

FgdGraphRecord::FgdGraphRecord(Int_t vecDist
            , Int_t disDiff
            , Double_t gradDiffAlowed)
        :  fvecDist(vecDist)
        , fdisDiff(disDiff)
        , fgradDiffAlowed(gradDiffAlowed)
{
    
}

FgdGraphRecord::~FgdGraphRecord()
{

}

FgdGraphRecord::FgdGraphRecord(const FgdGraphRecord& c)
        : fvecDist(c.fvecDist)
        , fdisDiff(c.fdisDiff)
        , fgradDiffAlowed(c.fgradDiffAlowed)
{
    
}

FgdGraphRecord& FgdGraphRecord::operator=(const FgdGraphRecord& c)
{
    fvecDist = c.fvecDist;
    fdisDiff = c.fdisDiff;
    fgradDiffAlowed = c.fgradDiffAlowed;
    return *this;
}

Bool_t FgdGraphRecord::operator==(const FgdGraphRecord& c)
{
    return (this->fvecDist == c.fvecDist
            && this->fdisDiff == c.fdisDiff);
}

} //superfgd
} //reconstruction
} //esbroot