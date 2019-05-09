 ################################################################################
 #    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    #
 #                                                                              #
 #              This software is distributed under the terms of the             # 
 #              GNU Lesser General Public Licence (LGPL) version 3,             #  
 #                  copied verbatim in the file "LICENSE"                       #
 ################################################################################
# Create a library called "libEsbEsb" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
${BASE_INCLUDE_DIRECTORIES}
#put here all directories where header files are located
${CMAKE_SOURCE_DIR}
#${CMAKE_SOURCE_DIR}/EsbGeometry
#${CMAKE_SOURCE_DIR}/EsbData
#${CMAKE_SOURCE_DIR}/esb
)

include_directories(${INCLUDE_DIRECTORIES})
include_directories(${SYSTEM_INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
${FAIRROOT_LIBRARY_DIR}
) 

link_directories( ${LINK_DIRECTORIES})

#Find all source files
file(GLOB sourceFiles ${CMAKE_SOURCE_DIR}/EsbGeometry/*.cxx)

set(SRCS
	${sourceFiles}
)

#Find all header files
file(GLOB headerFiles ${CMAKE_SOURCE_DIR}/EsbGeometry/*.h)

#GeometryLinkDef.h should be the last item in the list
list(REMOVE_ITEM headerFiles ${CMAKE_SOURCE_DIR}/EsbGeometry/GeometryLinkDef.h)

Set(HEADERS 
	${headerFiles}
)

Set(LINKDEF ${CMAKE_SOURCE_DIR}/EsbGeometry/GeometryLinkDef.h)
Set(LIBRARY_NAME EsbGeometry)
Set(DEPENDENCIES
    Base EsbData
)

GENERATE_LIBRARY()
