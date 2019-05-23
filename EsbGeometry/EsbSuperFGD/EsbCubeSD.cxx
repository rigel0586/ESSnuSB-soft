#include "EsbGeometry/EsbSuperFGD/EsbCubeSD.h"

#include "EsbGeometry/EsbSuperFGD/EsbFgdRunManager.h"
#include "EsbGeometry/EsbSuperFGD/EsbDetectorConstruction.h"

#include "FairLogger.h"

namespace esbroot {
namespace geometry {
namespace superfgd {

CubeSD::CubeSD(G4String name, G4String detector_sd_type)
    :G4VSensitiveDetector(name),fverbose(false), fdetector_sd_type(detector_sd_type) 
{
    fhitBuffer =  make_shared<data::superfgd::detector::CubeHit>();
}

CubeSD::~CubeSD()
{
}

void  CubeSD::Initialize(G4HCofThisEvent* HCE)
{

}

G4bool CubeSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ROHist)
{
    return ProcessTotalHits(aStep, ROHist);
}

/*
Process total hits sums up all the deposited energy and
trackingLenght before the information is exported.
This is done to be more memory efficient for large simulations.
*/
G4bool CubeSD::ProcessTotalHits(G4Step* aStep,G4TouchableHistory* ROHist)
{
    if(fverbose)
        LOG(debug) << " CubeSD::ProcessHitsSuperFGD ";

    //==========================================================================
    //              Get Runmanager and detector to analyze the step data
    //==========================================================================
    FgdRunManager* man = dynamic_cast<FgdRunManager*>(G4RunManager::GetRunManager());
    G4VUserDetectorConstruction* g4dc = const_cast<G4VUserDetectorConstruction*>(man->GetUserDetectorConstruction());
    FgdDetectorConstruction* detector = dynamic_cast<FgdDetectorConstruction*>(g4dc);
    if(detector==nullptr) 
        return false;

    bool ignoreEdep = (man!=nullptr && man->GetIgnoreEdep());

    //==============================================
    //              Get Step data
    //==============================================
    G4double edep = aStep->GetTotalEnergyDeposit();
    if (!ignoreEdep && edep == 0.)
        return false;

    //
    G4double nonIon = aStep->GetNonIonizingEnergyDeposit();
    G4ThreeVector hitPosition   = aStep->GetPreStepPoint()->GetPosition();
    G4ThreeVector postHitPosition = aStep->GetPostStepPoint()->GetPosition();
    G4double      kineticEnergy = aStep->GetTrack()->GetKineticEnergy();
    G4ThreeVector momentum = aStep->GetPreStepPoint()->GetMomentum();

    double time = aStep->GetPreStepPoint()->GetGlobalTime();
    double postTime = aStep->GetPostStepPoint()->GetGlobalTime();

    G4Track* track = aStep->GetTrack();
    int trackID = track->GetTrackID();
    int parentID = track->GetParentID();
    int pdg = track->GetDefinition()->GetPDGEncoding();
    G4double charge = track->GetDefinition()->GetPDGCharge();

    double trackLength = aStep->GetStepLength();
    
    if(man!=nullptr)
    {
        fwriter = man->GetFileWriter();
    }

    //================================================================
    //              Write data to be exported to file on end event
    //================================================================
    if(fwriter!=nullptr)
    {
        if(fverbose)
            LOG(debug) <<" ======  CubeSD:fverbose =====";

        fhitBuffer->SetEdep(edep);
        fhitBuffer->SetNonIonizingEnergyDeposit(nonIon);

        fhitBuffer->SetHitPostion(hitPosition.x(), hitPosition.y(), hitPosition.z());
        fhitBuffer->SetPostHitPostion(postHitPosition.x(), postHitPosition.y(), postHitPosition.z());
        fhitBuffer->SetKinEnergy(kineticEnergy);
        fhitBuffer->SetHitMomentum(momentum.x(), momentum.y(), momentum.z());

        fhitBuffer->SetTime(time);
        fhitBuffer->SetPostTime(postTime);
        fhitBuffer->SetTracklength(trackLength);

        fhitBuffer->SetTrackId(trackID);
        fhitBuffer->SetParentId(parentID);
        fhitBuffer->SetPdg(pdg);
        fhitBuffer->SetCharge(charge);

        //================================================================
        //     Check in which cube is the step position and accumulate it
        //================================================================
        data::superfgd::detector::FgdDetectorParameters& dp = detector->GetDetectorParams();
        fwriter->SumStep(*fhitBuffer.get(), dp);
    }
    else if(fverbose)
    {
        LOG(debug) << " FileWriter is nullptr ";
    }
    
    //================================================================
    //              Output verbose data
    //================================================================
    if(fverbose)
    {
        // Get Material 
        G4String thisVolume = aStep->GetTrack()->GetVolume()->GetName() ;
        G4String particleName = aStep->GetTrack()->GetDefinition()->GetParticleName();

        LOG(debug) << " ******************************* ";
        LOG(debug) << " HIT ";
        LOG(debug) << " Total_Energy_Deposit " << edep/eV;
        LOG(debug) << " Total_Energy_NonIoniziong " << nonIon/eV;
        LOG(debug) << " Pdg Code " << pdg;
        LOG(debug) << " Particle: " << particleName;
        LOG(debug) << " Volume: " << thisVolume;
        LOG(debug) << " TrackLength: " << trackLength / CLHEP::mm;
        LOG(debug) << " Charge: " << track->GetDefinition()->GetPDGCharge();
        LOG(debug) << " Energy (eV): " << kineticEnergy/eV;
        LOG(debug) << " POSITION (mm): " << hitPosition.x()/mm << hitPosition.y()/mm << hitPosition.z()/mm;
        LOG(debug) << " ******************************* ";
    }

    return true;
}

}   //superfgd
}   //geometry
}   //esbroot