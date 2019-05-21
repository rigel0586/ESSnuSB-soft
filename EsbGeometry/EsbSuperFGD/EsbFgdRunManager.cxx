#include "EsbGeometry/EsbSuperFGD/EsbFgdRunManager.h"

namespace esbroot {
namespace geometry {
namespace superfgd {

FgdRunManager::FgdRunManager():G4RunManager(),ffileWriter(nullptr),fignoreEdep(false)
{
}

FgdRunManager::~FgdRunManager()
{
}

}   //superfgd
}   //geometry
}   //esbroot