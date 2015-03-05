
cmake_minimum_required(VERSION 2.8)

PROJECT(readEGS_in_VTK_proj)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

#you need to set this to the root directory of our attila source directory
SET(ATTILA_SRC_DIR /home/jkmikell/src/meshGenGitHub)


##BOOST SPECIFIC STUFF########################################################
#you can tell cmake where to find boost here if installed in nonstandard place
SET(BOOST_ROOT "~/apps/boost_1_53_0")
SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "${BOOST_ROOT}/include")
SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "${BOOST_ROOT}/lib")
#boost settings/find boost
SET(Boost_USE_STATIC_LIBS OFF)
SET(Boost_USE_MULTITHREADED ON)
SET(Boost_USE_STATIC_RUNTIME OFF)
FIND_PACKAGE(Boost 1.53 COMPONENTS program_options filesystem system regex unit_test_framework)
IF (Boost_FOUND)
        MESSAGE("Boost found!!!")
        INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
        LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})
        ADD_DEFINITIONS( "-DHAS_BOOST" )
ENDIF()
###END BOOST SPECIFIC STUFF#################################################

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
#############################################################################
#=========================================================

# Locate sources and headers for this project

# vtk read/write mesh fileformat test
SET(MYVTKDIR ${ATTILA_SRC_DIR}/src/vtkMeshFormat)

FILE(GLOB headers
        ${ATTILA_SRC_DIR}/include/*.hpp)

#Locate Attila files
FILE(GLOB attilaSources
        ${ATTILA_SRC_DIR}/src/Attila*.cpp)

FILE(GLOB attilaHeaders
        ${ATTILA_SRC_DIR}/include/Attila*.hpp)

FILE(GLOB testSources
        ${TEST_SOURCE_DIR}/*_test.cpp)

MESSAGE("")
MESSAGE("attila headers: ")
foreach(TEMPVAR ${attilaHeaders})
        MESSAGE("${TEMPVAR}")
endforeach()
MESSAGE("")

MESSAGE("boost library directories:")
foreach(TEMPVAR ${Boost_LIBRARY_DIRS})
        MESSAGE("${TEMPVAR}")
endforeach()
MESSAGE("")

MESSAGE("boost include directory:")
foreach(TEMPVAR ${Boost_INCLUDE_DIR})
        MESSAGE("${TEMPVAR}")
endforeach()
MESSAGE("")

MESSAGE("boost libraries:")
foreach(TEMPVAR ${Boost_LIBRARIES})
        MESSAGE("${TEMPVAR}")
endforeach()
MESSAGE("")

###########################################################
#############################################################################
# Setup include directory for this project
INCLUDE_DIRECTORIES(
        "${ATTILA_SRC_DIR}/include"
        "${Boost_INCLUDE_DIRS}")
#############################################################################

SET(LIB_VAR egsphantff)
ADD_LIBRARY("${LIB_VAR}Lib" ${ATTILA_SRC_DIR}/src/${LIB_VAR}.cpp)
TARGET_LINK_LIBRARIES("${LIB_VAR}Lib" xyzBoundArrLib ${Boost_LIBRARIES})

SET(LIB_VAR xyzBoundArr)
ADD_LIBRARY("${LIB_VAR}Lib" ${ATTILA_SRC_DIR}/src/${LIB_VAR}.cpp)        
TARGET_LINK_LIBRARIES("${LIB_VAR}Lib" AttilaUtilitiesLib ${Boost_LIBRARIES})

SET(LIB_VAR AttilaUtilities)
add_library("${LIB_VAR}Lib" ${ATTILA_SRC_DIR}/src/${LIB_VAR}.cpp)
target_link_libraries(${LIB_VAR}Lib ${Boost_LIBRARIES})
SET(LIB_VAR "")
##############################################################################

SET(LIB_VAR node_selection)
ADD_LIBRARY("${LIB_VAR}Lib" ${LIB_VAR}.cpp)
TARGET_LINK_LIBRARIES("${LIB_VAR}Lib" ${Boost_LIBRARIES})

add_executable(readEGS_in_VTK MACOSX_BUNDLE readEGS_in_VTK)
add_executable(generate_nodes_from_prob_vertices generate_nodes_from_prob_vertices)


if(VTK_LIBRARIES)
  target_link_libraries(readEGS_in_VTK ${VTK_LIBRARIES} egsphantffLib)
else()
  target_link_libraries(readEGS_in_VTK vtkHybrid vtkWidgets egsphantffLib)
endif()
