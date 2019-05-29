#ifndef ESBROOT_ESBGEOMETRY_SUPERFGD_RUN_ACTION_H
#define ESBROOT_ESBGEOMETRY_SUPERFGD_RUN_ACTION_H 1

#include <G4UserRunAction.hh>
#include <globals.hh>

class G4Run;

namespace esbroot {
namespace geometry {
namespace superfgd {

class FgdRunAction : public G4UserRunAction
{
public:

    /** Default constructor **/
    FgdRunAction();

    /** Destructor **/
    virtual ~FgdRunAction();

    /** Base method to begin the run action
     *@param G4Run - current run
     **/
    virtual void BeginOfRunAction(const G4Run*);

    /** Base method to end the run action
     *@param G4Run - current run
     **/
    virtual void   EndOfRunAction(const G4Run*); 
};

}   //superfgd
}   //geometry
}   //esbroot

#endif