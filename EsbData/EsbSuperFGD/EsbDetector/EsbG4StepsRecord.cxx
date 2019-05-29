#include "EsbData/EsbSuperFGD/EsbDetector/EsbG4StepsRecord.h"

namespace esbroot {

namespace data {

namespace superfgd {

namespace detector {

G4StepsRecord::G4StepsRecord()
    : fbin_x(-1),fbin_y(-1),fbin_z(-1), fkey(-1)
    , fstartTime(0.), fendTime(0.), fedep(0.)
    , ftrackLenght(0.),fnonIon(0.)
{
}

G4StepsRecord::~G4StepsRecord()
{
}


bool G4StepsRecord::operator==(const G4StepsRecord& rhs) const
{
    return rhs.GetBinX()==this->GetBinX()
            && rhs.GetBinY()==this->GetBinY()
            && rhs.GetBinZ()==this->GetBinZ();
}
   

void G4StepsRecord::Copy(TObject& object) const
{
    G4StepsRecord* step = dynamic_cast<G4StepsRecord*>(&object);

    step->SetBinPosition(this->fbin_x,this->fbin_y,this->fbin_z);
    step->SetStartTime(this->fstartTime);
    step->SetEndTime(this->fendTime);
    step->SetEdep(this->fedep);
    step->SetTrackLenght(this->ftrackLenght);
    step->SetNonIon(this->fnonIon);
    step->SetTrackIds(const_cast<std::set<Int_t>&>(this->ftrackIds));
    step->SetParentIds(const_cast<std::set<Int_t>&>(this->fparentIds));
    step->SetPdgs(const_cast<std::set<Int_t>&>(this->fpdgs));
}

}   // detector
}   // superfgd
}   // data
}   // esbroot