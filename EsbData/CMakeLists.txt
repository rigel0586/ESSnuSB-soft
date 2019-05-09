 # Create a library called "libEsbData" 

set(INCLUDE_DIRECTORIES
${BASE_INCLUDE_DIRECTORIES}
#${CMAKE_SOURCE_DIR}/EsbData
${CMAKE_SOURCE_DIR}
)

include_directories( ${INCLUDE_DIRECTORIES})
include_directories(${SYSTEM_INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
${FAIRROOT_LIBRARY_DIR}
#${CMAKE_SOURCE_DIR}/EsbData
)
 
link_directories( ${LINK_DIRECTORIES})

#Find all source files
file(GLOB sourceFiles ${CMAKE_SOURCE_DIR}/EsbData/*.cxx)

set(SRCS
    ${sourceFiles}
)

#Find all header files
file(GLOB headerFiles ${CMAKE_SOURCE_DIR}/EsbData/*.h)

#DataLinkDef.h should be the last item in the list
list(REMOVE_ITEM headerFiles ${CMAKE_SOURCE_DIR}/EsbData/DataLinkDef.h)

Set(HEADERS 
    ${headerFiles}
)

Set(LINKDEF ${CMAKE_SOURCE_DIR}/EsbData/DataLinkDef.h)
Set(LIBRARY_NAME EsbData)
Set(DEPENDENCIES Base EG Physics Core)

GENERATE_LIBRARY()
