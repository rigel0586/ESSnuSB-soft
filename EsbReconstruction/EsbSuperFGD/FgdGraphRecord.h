#ifndef ESBROOT_ESBDRECONSTRUCTION_FGD_GRAPH_RECORD_H
#define ESBROOT_ESBDRECONSTRUCTION_FGD_GRAPH_RECORD_H

#include "TObject.h"
#include <TVector3.h>
#include <vector>
#include <algorithm>

namespace esbroot {
namespace reconstruction {
namespace superfgd {

class FgdGraphRecord : public TObject
{

public:

    /** Default constructor **/  
    FgdGraphRecord();

    /** Constructor with argument
    *@param vecDist  - number of cubes as a track segment to calculate the direction vector
    *@param disDiff  - number of distance cubes to calculate the difference vector
    *@param gradDiffAlowed  - allowable difference in degrees
    **/  
    FgdGraphRecord(Int_t vecDist
            , Int_t disDiff
            , Double_t gradDiffAlowed);

    ~FgdGraphRecord();

    FgdGraphRecord(const FgdGraphRecord& c);

    FgdGraphRecord& operator=(const FgdGraphRecord& c);

    Bool_t operator==(const FgdGraphRecord& c);

    TVector3 operator-(const FgdGraphRecord& c);


    Int_t fvecDist;
    Int_t fdisDiff;
    Double_t fgradDiffAlowed;

private:

    ClassDef(FgdGraphRecord, 2);
};

} //superfgd
} //reconstruction
} //esbroot

#endif