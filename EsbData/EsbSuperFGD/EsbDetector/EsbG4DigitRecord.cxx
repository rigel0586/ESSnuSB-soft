#include "EsbData/EsbSuperFGD/EsbDetector/EsbG4DigitRecord.h"

ClassImp(esbroot::data::superfgd::detector::G4DigitRecord)

namespace esbroot {

namespace data {

namespace superfgd {

namespace detector {
    
ostream & operator << (ostream & stream, const G4DigitRecord & rec) 
{
    return stream;
}

G4DigitRecord::G4DigitRecord()
:TObject()
{
    freadouts = new TClonesArray("EsbDigitReadout");
    freadouts->SetOwner(true);
    Init();
}


G4DigitRecord::~G4DigitRecord() 
{
    Clear("C");
    delete freadouts;
}

DigitReadout* G4DigitRecord::GetDigitReadout(Int_t index) const
{
    if ((index >= 0) && (index < freadouts->GetEntriesFast())) 
    {
        DigitReadout* ro = (DigitReadout*)(*freadouts)[index];
        if (ro && ro!=nullptr && ro!=NULL)
        {
            return ro;
        }
    }
    
    return NULL;
}


unsigned int G4DigitRecord::NumDigitReadouts() const 
{
    return freadouts->GetEntries();
}

void G4DigitRecord::AddDigitReadout(const DigitReadout& ro)
{
    unsigned int pos = freadouts->GetEntries();
    new((*freadouts)[pos]) DigitReadout(ro);
}


void G4DigitRecord::Copy(TObject& object) const 
{
    G4DigitRecord* record = dynamic_cast<G4DigitRecord*>(&object);
    record->Reset();

    // copy cubeHits
    TClonesArray& srcRo = *(this->freadouts);;
    TClonesArray& destRo = *(record->freadouts);

    int nc = freadouts->GetEntries();
    for (int i = 0; i < nc; i++) 
    {
        new(destRo[i]) DigitReadout(*((DigitReadout*)srcRo[i]));
    }
}

void G4DigitRecord::Clear(Option_t* /*option*/) 
{
    freadouts->Clear("C");
}

void G4DigitRecord::Print(Option_t* option) const 
{
    
}

void G4DigitRecord::Init() 
{
}

void G4DigitRecord::Reset() 
{
    Clear();
    Init();
}

}   // detector
}   // superfgd
}   // data
}   // esbroot