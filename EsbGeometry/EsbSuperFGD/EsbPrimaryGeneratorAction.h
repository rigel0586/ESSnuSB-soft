#ifndef ESBROOT_ESBGEOMETRY_SUPERFGD_PRIMARY_GENERATOR_ACTION_H
#define ESBROOT_ESBGEOMETRY_SUPERFGD_PRIMARY_GENERATOR_ACTION_H 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

#include "G4ParticleDefinition.hh"
#include "Framework/GHEP/GHepParticle.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "Framework/Ntuple/NtpMCEventRecord.h"
#include <TFile.h>
#include <TChain.h>

class G4Event;

namespace esbroot {
namespace geometry {
namespace superfgd {

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:

    /** Constructor
     *@param eventsFile - file containing the generated events to read from
     *@param phys - physical list of the allowed events and particle types
     **/
    PrimaryGeneratorAction(std::string eventsFile);

    /** Destructor **/
    virtual ~PrimaryGeneratorAction();

    /** method from the base class **/
    virtual void GeneratePrimaries(G4Event*);  
    
    /** Set verbose mode for the primary generator
     *@param mode - mode to log
     **/
    void SetVerboseMode(bool mode){fverboseMode = mode;}

    /** Returns number of events in the eventsFile **/
    int GetNumberOfEvents() {return fchainNEvents;}

    /** Returns true if there are more events to read, false otherwise
     * counter is incremented on each call to ReadEvent **/
    bool HasMoreEvents() const {return (fnextEventNo < fchainNEvents);}
  
  private:
    
    /** Initializes the TChain from the eventsFile **/
    void InitializeTChain();

    /** Adds a particle to the initial vertex of the interaction **/
    void AddParticleToVertex(G4PrimaryVertex* v, genie::GHepParticle *p);

    /** Returns an event from the top of the stack **/
    const genie::NtpMCEventRecord* PopEvent();

    /** Read the next event from the eventsFile **/
    bool ReadEvent();

    /** Contains a pointer to the genie generated event **/
    genie::NtpMCEventRecord* fbuffer;

    /** Copy of the fbuffer, to be return (if modified) from Pop **/
    genie::NtpMCEventRecord* foldbuffer;

    /** Chain of events in the file **/
    TChain* fchain;

    /** verbose mode **/
    bool fverboseMode;

    /** Total number of events **/
    int fchainNEvents;

    /** Counter for the next event **/
    int fnextEventNo;

    /** full path of the file containing the events **/
    std::string feventsFile;
};

}   //superfgd
}   //geometry
}   //esbroot

#endif