 # Create a library called "libEsbData" 

message(STATUS "Using Genie : $ENV{GENIE}")
message(STATUS "Using CLHEP_INCLUDE_DIR : ${CLHEP_INCLUDE_DIR}")

set(INCLUDE_DIRECTORIES
    ${BASE_INCLUDE_DIRECTORIES}
    #${CMAKE_SOURCE_DIR}/EsbData
    ${CMAKE_SOURCE_DIR}
    ${CLHEP_INCLUDE_DIR}
    $ENV{GENIE}/src
)

include_directories( ${INCLUDE_DIRECTORIES})
include_directories(${SYSTEM_INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
    ${ROOT_LIBRARY_DIR}
    ${FAIRROOT_LIBRARY_DIR}
    #${CMAKE_SOURCE_DIR}/EsbData
    ${CLHEP_LIB_DIR}
    $ENV{GENIE}/lib
)
 
link_directories( ${LINK_DIRECTORIES})

file(GLOB fgdDetectorSources ${CMAKE_SOURCE_DIR}/EsbData/EsbSuperFGD/EsbDetector/*.cxx)

set(SRCS
    FgdDetectorPoint.cxx
    MCTrack.cxx
    Stack.cxx
    PMTubeHit.cxx
    WCDetectorPoint.cxx
    ${fgdDetectorSources}
)

file(GLOB fgdDetectorHeaders ${CMAKE_SOURCE_DIR}/EsbData/EsbSuperFGD/EsbDetector/*.h)

Set(HEADERS 
    ${CMAKE_SOURCE_DIR}/EsbData/FgdDetectorPoint.h
    ${CMAKE_SOURCE_DIR}/EsbData/MCTrack.h
    ${CMAKE_SOURCE_DIR}/EsbData/PMTubeHit.h
    ${CMAKE_SOURCE_DIR}/EsbData/Stack.h
    ${CMAKE_SOURCE_DIR}/EsbData/WCDetectorPoint.h
    ${fgdDetectorHeaders}
)

Set(LINKDEF ${CMAKE_SOURCE_DIR}/EsbData/DataLinkDef.h)
Set(LIBRARY_NAME EsbData)
Set(DEPENDENCIES Base EG Physics Core)

GENERATE_LIBRARY()