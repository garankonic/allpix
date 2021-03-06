# Locate Pythia8 library
# in a directory defined via PYTHIA8 environment variable
#
# Defines:
#  PYTHIA8_FOUND
#  PYTHIA8_INCLUDE_DIR
#  PYTHIA8_LIBRARY
#  PYTHIA8_HEPMC_LIBRARY
#  PYTHIA8_HAPDF_LIBRARY

find_path(PYTHIA8_INCLUDE_DIR Pythia8/Pythia.h Pythia8/Pythia8ToHepMC.h
          HINTS $ENV{PYTHIA8LOCATION}/include ${PYTHIA8LOCATION}/include
          $ENV{PYTHIA8}/include ${PYTHIA8}/include)

# message(STATUS PYTHIA8_INCLUDE_DIR ${PYTHIA8_INCLUDE_DIR} )

set(ENV{PYTHIA8} "/home/maren/GEANT4/pythia8185_install")

find_library(PYTHIA8_LIBRARY NAMES pythia8
			     HINTS $ENV{PYTHIA8LOCATION}/lib ${PYTHIA8LOCATION}/lib
			     HINTS $ENV{PYTHIA8}/lib/archive ${PYTHIA8}/lib/archive)

# message(STATUS PYTHIA8_LIBRARY ${PYTHIA8_LIBRARY} )

find_library(PYTHIA8_HEPMC_LIBRARY NAMES pythia8tohepmc
			     HINTS $ENV{PYTHIA8LOCATION}/lib ${PYTHIA8LOCATION}/lib
			     HINTS $ENV{PYTHIA8}/lib/archive ${PYTHIA8}/lib/archive)

# message(STATUS PYTHIA8_HEPMC_LIBRARY ${PYTHIA8_HEPMC_LIBRARY} )

find_library(PYTHIA8_HAPDF_LIBRARY NAMES lhapdfdummy
			     HINTS $ENV{PYTHIA8LOCATION}/lib ${PYTHIA8LOCATION}/lib
			     HINTS $ENV{PYTHIA8}/lib/archive ${PYTHIA8}/lib/archive)

# message(STATUS PYTHIA8_HAPDF_LIBRARY ${PYTHIA8_HAPDF_LIBRARY} )

# handle the QUIETLY and REQUIRED arguments and set PYTHIA8_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Pythia8 DEFAULT_MSG PYTHIA8_LIBRARY)

mark_as_advanced(PYTHIA8_FOUND PYTHIA8_LIBRARY PYTHIA8_HEPMC_LIBRARY PYTHIA8_HAPDF_LIBRARY)
