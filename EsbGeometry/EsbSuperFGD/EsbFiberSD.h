#ifndef ESBROOT_ESBGEOMETRY_SUPERFGD_FIBER_SD_H
#define ESBROOT_ESBGEOMETRY_SUPERFGD_FIBER_SD_H 1

#include <memory>

#include <G4VSensitiveDetector.hh>

#include "EsbData/EsbSuperFGD/EsbDetector/EsbFiberHit.h"
#include "EsbGeometry/EsbSuperFGD/EsbFileWriter.h"

namespace esbroot {
namespace geometry {
namespace superfgd {

class FiberSD : public G4VSensitiveDetector
{
public:
    /** Default constructor **/
    FiberSD(G4String name);

    /** Destructor **/
    virtual ~FiberSD();

    /** Process the hit from the geant4 simulation in the fiber
     * extacting simulation data
     *@param astep - step in the simulation
     *@param ROHist - geant4 history
     **/
    G4bool ProcessHits(G4Step* astep,G4TouchableHistory* ROHist);

    /** Initialization of Fiber Sensitive Detector
     *@param HCE - configuration parameter (not used)
     **/
    void Initialize(G4HCofThisEvent* HCE);

    /** Set verbosity level
     *@param verbose - true to verbose
     **/
    void SetVerbose(bool verbose) {fverbose = verbose;}

private:
    std::shared_ptr<data::superfgd::detector::FiberHit> fhitBuffer;
    std::shared_ptr<FileWriter> fwriter;
    bool fverbose;
};

}   //superfgd
}   //geometry
}   //esbroot

#endif