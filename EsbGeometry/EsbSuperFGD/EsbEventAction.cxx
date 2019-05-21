#include "EsbGeometry/EsbSuperFGD/EsbEventAction.h"
#include "EsbGeometry/EsbSuperFGD/EsbFileWriter.h"
#include "EsbGeometry/EsbSuperFGD/EsbFgdRunManager.h"

#include "EsbGeometry/EsbSuperFGD/EsbPrimaryGeneratorAction.h"
#include "EsbGeometry/EsbSuperFGD/EsbParticleGunPrimaryGenerator.h"

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4RunManager.hh"


namespace esbroot {
namespace geometry {
namespace superfgd {

FgdEventAction::FgdEventAction() 
{
}
 
FgdEventAction::~FgdEventAction() 
{
}
 
void FgdEventAction::BeginOfEventAction(const G4Event* event) 
{
}
 
void FgdEventAction::EndOfEventAction(const G4Event* event)
{
    G4RunManager *runmanager = G4RunManager::GetRunManager();
    G4VUserPrimaryGeneratorAction* prGA = const_cast<G4VUserPrimaryGeneratorAction*>(runmanager->GetUserPrimaryGeneratorAction());
    
    PrimaryGeneratorAction* nfPrga = dynamic_cast<PrimaryGeneratorAction*>(prGA);
    ParticleGunPrimaryGenerator* nfPrGga = dynamic_cast<ParticleGunPrimaryGenerator*>(prGA);

    std::shared_ptr<FileWriter> writer = ((FgdRunManager*)G4RunManager::GetRunManager())->GetFileWriter();
    writer->WriteHit();

    if(nfPrga!=nullptr && !nfPrga->HasMoreEvents()
        || nfPrGga!=nullptr && !nfPrGga->HasMoreEvents())
    {
        G4cout << "FgdEventAction: No more events, aborting run ..." << G4endl;
        runmanager->AbortRun();
    }
}

}   //superfgd
}   //geometry
}   //esbroot