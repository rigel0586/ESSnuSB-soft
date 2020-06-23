#include "EsbReconstruction/EsbSuperFGD/FgdMCEventRecord.h"

#include "Framework/ParticleData/PDGCodes.h"
#include "FairLogger.h"

#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <memory>
#include <math.h>
#include <bits/stdc++.h>

namespace esbroot {
namespace reconstruction {
namespace superfgd {


// -----   Constructors and Destructors  -----------------------------------

FgdMCEventRecord::FgdMCEventRecord(std::string eventData) 
    : TObject(), feventData(eventData), fnuPdg(0)
        , fNuEnergy(0.), fvertex(TVector3(0,0,0)), fIsWeakCC(false)
        , fIsWeakNC(false), fIsQuasiElastic(false), fIsDeepInElastic(false)
        , fIsResonant(false), fIsCoherentProduction(false), fIsMEC(false)
        , fIsNuElectronElastic(false), fIsElectronScattering(false)
        , fIsPrimaryMuon(false)
        , fIsPrimaryElectron(false), fPrimaryMuonMom(TVector3(0,0,0))
        , fPrimaryElectronMom(TVector3(0,0,0)), fMuonTrackLength(0.)
        , fIsMuonExiting(false), fMuonExitMomentum(TVector3(0,0,0))
        , fMuonPolarAngle(0.), fMuonAzumAngle(0.), fElectronNumOfExitingParticles(0)
        , fGenfitMom(TVector3(0,0,0)), fMC_GentFitError(0.), fProtonEdep(0.)
        , fHasHits(false), fDeDx(0.), fEdepDeDx(0.)
{
    Init();
}

FgdMCEventRecord::FgdMCEventRecord() 
    : feventData("")
{
}

FgdMCEventRecord::FgdMCEventRecord(const FgdMCEventRecord& c)
{
    this->feventData = c.feventData;
    
    this->fMuonTrackLength = c.fMuonTrackLength;
    this->fIsMuonExiting = c.fIsMuonExiting;
    this->fMuonExitMomentum = c.fMuonExitMomentum;
    this->fMuonPolarAngle = c.fMuonPolarAngle;
    this->fMuonAzumAngle = c.fMuonAzumAngle;
    this->fHadronEdep = c.fHadronEdep;
    this->fElectronNumOfExitingParticles = c.fElectronNumOfExitingParticles;
    this->fGenfitMom = c.fGenfitMom;
    this->fMC_GentFitError = c.fMC_GentFitError;
    this->fProtonEdep = c.fProtonEdep;
    this->fHasHits = c.fHasHits;
    this->fDeDx = c.fDeDx;
    this->fEdepDeDx = c.fEdepDeDx;

    this->Init();
}

FgdMCEventRecord::~FgdMCEventRecord()
{
}
// -------------------------------------------------------------------------

// -----   Public methods   ------------------------------------------------

void FgdMCEventRecord::SetEventData(std::string data)
{
    feventData = data;

    Init();
}

Int_t FgdMCEventRecord::GetNuPdg(void)
{
    return fnuPdg;
}

Double_t FgdMCEventRecord::GetNuE(void)
{
    return fNuEnergy;
}

TVector3 FgdMCEventRecord::GetVertex(void)
{
    return fvertex;
}

Bool_t FgdMCEventRecord::IsWeakCC(void)
{
    return fIsWeakCC;
}

Bool_t FgdMCEventRecord::IsWeakNC(void)
{
    return fIsWeakNC;
}

Bool_t FgdMCEventRecord::IsQuasiElastic(void)
{
    return fIsQuasiElastic;
}

Bool_t FgdMCEventRecord::IsDeepInElastic(void)
{
    return fIsDeepInElastic;
}

Bool_t FgdMCEventRecord::IsResonant(void)
{
    return fIsResonant;
}

Bool_t FgdMCEventRecord::IsCoherentProduction(void)
{
    return fIsCoherentProduction;
}

Bool_t FgdMCEventRecord::IsMEC(void)
{
    return fIsMEC;
}

Bool_t FgdMCEventRecord::IsNuElectronElastic(void)
{
    return fIsNuElectronElastic;
}  

Bool_t FgdMCEventRecord::IsElectronScattering(void)
{
    return fIsElectronScattering;
}

Bool_t FgdMCEventRecord::HasHits(void)
{
    return fHasHits;
}


Bool_t FgdMCEventRecord::IsPrimaryLeptonMuon()
{
    return fIsPrimaryMuon;
}

TVector3 FgdMCEventRecord::GetMuonMom()
{
    return fPrimaryMuonMom;
}

Double_t FgdMCEventRecord::GetMuonTrackLength()
{
    return fMuonTrackLength;
}


Double_t FgdMCEventRecord::GetMuonTrackLengthOrigin()
{
    return fMuonTrackLengthOrigin;
}

Bool_t FgdMCEventRecord::IsPrimaryLeptonElectron()
{
    return fIsPrimaryElectron;
}

TVector3 FgdMCEventRecord::GetElectronMom()
{
    return fPrimaryElectronMom;
}

Double_t FgdMCEventRecord::GetMuonPolarAngle()
{
    return fMuonPolarAngle;
}

Double_t FgdMCEventRecord::GetMuonAzumuteAngle()
{
    return fMuonAzumAngle;
}

void FgdMCEventRecord::SetMuonExitMom(TVector3 exitMom)
{
    fMuonExitMomentum = exitMom;

    Double_t radToAngle = 180/TMath::Pi();
    fMuonPolarAngle = fMuonExitMomentum.Theta() * radToAngle;
    fMuonAzumAngle = fMuonExitMomentum.Phi() * radToAngle;
}


FgdMCEventRecord& FgdMCEventRecord::operator=(const FgdMCEventRecord& c)
{
    this->feventData = c.feventData;

    this->fMuonTrackLength = c.fMuonTrackLength;
    this->fIsMuonExiting = c.fIsMuonExiting;
    this->fMuonExitMomentum = c.fMuonExitMomentum;
    this->fMuonPolarAngle = c.fMuonPolarAngle;
    this->fMuonAzumAngle = c.fMuonAzumAngle;
    this->fHadronEdep = c.fHadronEdep;
    this->fElectronNumOfExitingParticles = c.fElectronNumOfExitingParticles;
    this->fGenfitMom = c.fGenfitMom;
    this->fMC_GentFitError = c.fMC_GentFitError;
    this->fProtonEdep = c.fProtonEdep;
    this->fHasHits = c.fHasHits;
    this->fDeDx = c.fDeDx;
    this->fEdepDeDx = c.fEdepDeDx;

    this->Init();
    return *this;
}

const std::vector<std::pair<Int_t, TVector3>>& FgdMCEventRecord::GetPrimaryParticles()
{
    Init();

    return fPrimaryParticles;
}

void FgdMCEventRecord::PrintData(std::ostream & stream)
{
    stream << feventData;
}

void FgdMCEventRecord::ReadEventData()
{
    this->Init();
}

// -------------------------------------------------------------------------


// -----   Protected methods   ---------------------------------------------
void FgdMCEventRecord::Init()
{
    fDataTokens.clear();

    try
    {        
        const char spaceChar(' ');

        std::istringstream ss(feventData);
        std::string token;
        while(std::getline(ss, token, spaceChar))
        {
            if(!token.empty())
            {
                fDataTokens.push_back(token);
            }
        }
        InitMembers();
    }
    catch(const std::exception& e)
    {
        LOG(fatal) << e.what();
    }
}

void FgdMCEventRecord::InitMembers()
{
    fnuPdg = std::stoi(fDataTokens[Data::NEUTRINO_PDG]);
    fNuEnergy = std::stod(fDataTokens[Data::NEUTRINO_ENERGY]);

    Double_t x = std::stod(fDataTokens[Data::VERTEX_POSTION_X]);
    Double_t y = std::stod(fDataTokens[Data::VERTEX_POSTION_Y]);
    Double_t z = std::stod(fDataTokens[Data::VERTEX_POSTION_Z]);

    fvertex = TVector3(x,y,z);

    Int_t isWeakCC = std::stoi(fDataTokens[Data::IS_WEAK_CC]);
    fIsWeakCC = (isWeakCC==1);

    Int_t isWeakNC = std::stoi(fDataTokens[Data::IS_WEACK_NC]);
    fIsWeakNC = (isWeakNC==1);

    Int_t isQuasiElastic = std::stoi(fDataTokens[Data::IS_QUASI_ELASTIC]);
    fIsQuasiElastic = (isQuasiElastic==1);

    Int_t isDIS = std::stoi(fDataTokens[Data::IS_DEEP_INELASRIC]);
    fIsDeepInElastic = (isDIS==1);

    Int_t isResonant = std::stoi(fDataTokens[Data::IS_RESONANT]);
    fIsResonant = (isResonant==1);

    Int_t isCohPro = std::stoi(fDataTokens[Data::IS_COHERENT_PRODUCTION]);
    fIsCoherentProduction = (isCohPro==1);

    Int_t isMec = std::stoi(fDataTokens[Data::IS_MEC]);
    fIsMEC = (isMec==1);

    Int_t isNuElEl = std::stoi(fDataTokens[Data::IS_NuElectron_Elastic]);
    fIsNuElectronElastic = (isNuElEl==1); 

    Int_t isElEl = std::stoi(fDataTokens[Data::IS_ELECTRON_SCATTERING]);
    fIsElectronScattering = (isElEl==1); 

    // Initalize primary particles
    size_t particleDataSize = 4;

    size_t particlePdg = 0;
    size_t particleMomX = 1;
    size_t particleMomY = 2;
    size_t particleMomZ = 3;

    fPrimaryParticles.clear();
    for(size_t i = Data::PRIMARY_PARTICLES; i < fDataTokens.size(); i+=particleDataSize)
    {
        Int_t pdg = std::stoi(fDataTokens[i + particlePdg]);
        Double_t momX = std::stod(fDataTokens[i + particleMomX]);
        Double_t momY = std::stod(fDataTokens[i + particleMomY]);
        Double_t momZ = std::stod(fDataTokens[i + particleMomZ]);

        fPrimaryParticles.push_back(
                            std::pair<Int_t, TVector3>( pdg,    TVector3(momX , momY, momZ) )
                            );
    }

    for(size_t i = 0; i< fPrimaryParticles.size(); ++i)
    {
        std::pair<Int_t, TVector3>& particle = fPrimaryParticles[i];
        if(particle.first == genie::kPdgMuon || particle.first == genie::kPdgAntiMuon)
        {
            fIsPrimaryMuon = true;
            fPrimaryMuonMom = particle.second;
            break;
        }

        if(particle.first == genie::kPdgElectron || particle.first == genie::kPdgPositron)
        {
            fIsPrimaryElectron = true;
            fPrimaryElectronMom = particle.second;
            break;
        }
    }
}

// -------------------------------------------------------------------------

}// namespace superfgd
}// namespace reconstruction
}// namespace esbroot