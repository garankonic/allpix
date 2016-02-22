#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(ALLPIX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})
# ROOT
include(FindROOT.cmake)

set(INCLUDE_DIRECTORIES
	${ROOT_INCLUDE_DIR} 
)
include_directories( ${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
	${ROOT_LIBRARY_DIR}
)
link_directories( ${LINK_DIRECTORIES} )

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Find HepMC (required package)
#
find_package(HepMC REQUIRED)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Find Pythia8 (required package)
#
find_package(Pythia8 REQUIRED)
if(PYTHIA8_FOUND)
  message(STATUS "G4 Examples: Pythia8 found. --> HepMCEx03 example with Pythia8 enabled.") 
  add_definitions(-DG4LIB_USE_PYTHIA8) 
else()
  message(STATUS "G4 Examples: Pythia8 -- NOT found.") 
  set(PYTHIA8_LIBRARY "")
  set(PYTHIA8_HEPMC_LIBRARY "")
  set(PYTHIA8_HAPDF_LIBRARY "")
endif()

#set(HEPMC_INCLUDE_DIR "/home/maren/GEANT4/HepMC-2.06.09_install/include")

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR}
		    ${HEPMC_INCLUDE_DIR}
                    ${PYTHIA8_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/SelDict.cc ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

ROOT_GENERATE_DICTIONARY(HepMC include/HepMCG4Classes.hh  LINKDEF include/HepMCG4LinkDef.h )
add_library(HepMCdict SHARED ${sources} HepMC.cxx)
target_link_libraries(HepMCdict ${HEPMC_LIBRARIES} ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(allpix allpix.cc ${sources} ${headers})
target_link_libraries(allpix ${Geant4_LIBRARIES} ${HEPMC_LIBRARIES}
                       ${PYTHIA8_LIBRARY}
                       ${PYTHIA8_HEPMC_LIBRARY} ${PYTHIA8_HAPDF_LIBRARY}
                       HepMCdict ${ROOT_LIBRARIES})

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B5. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#

file(GLOB MACROLIST RELATIVE ${PROJECT_SOURCE_DIR}  "macros/*.in")
file(GLOB GDMLLIST RELATIVE ${PROJECT_SOURCE_DIR}  "models/*.gdml")

set(ALLPIX_SCRIPTS
    ${MACROLIST}
    ${GDMLLIST}
    macros/oneDetector_vis.in
    macros/mpx3rx_vis.in
    macros/mpx3rx.in
	models/pixeldetector.xml
	models/clicpix_box_materials.xml
	share/Erf.dat
	share/GDML_EUDETAlHolder/EUDETAlHolder.gdml
	share/GDML_EUDETAlHolder/EUDETAlHolder_materials.xml
  )

foreach(_script ${ALLPIX_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# rootcling - dictionaries
#add_custom_target(TARGET rootcling PRE_BUILD 
#COMMAND	@echo "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")

add_custom_target(rootcling_build ALL DEPENDS ALLPIX)

add_custom_command(TARGET rootcling_build PRE_BUILD
WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
COMMAND rootcling -v0 -f SelDict.cc -c -p -I./include
AllPix_Hits_WriteToEntuple.h 
AllPix_Frames_WriteToEntuple.h allpix_dm.h 
AllPixDigitAnimation.hh LinkDef.h
COMMAND mv SelDict.cc ./src/SelDict.cc
COMMAND mv SelDict_rdict.pcm ${PROJECT_BINARY_DIR}/SelDict_rdict.pcm)



#rootcint -v0 -f SelDict.cc -c -p -I./include 
#		AllPix_Hits_WriteToEntuple.h 
#		AllPix_Frames_WriteToEntuple.h allpix_dm.h 
#		AllPixDigitAnimation.hh LinkDef.h)

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(ALLPIX DEPENDS allpix)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS allpix DESTINATION bin)

