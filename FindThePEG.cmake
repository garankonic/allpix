# - Locate ThePEG library
# in a directory defined via  THEPEG_ROOT_DIR or THEPEG_DIR environment variable
# Defines:
#
#  THEPEG_FOUND
#  THEPEG_INCLUDE_DIR
#  THEPEG_INCLUDE_DIRS (not cached)
#  THEPEG_LIBRARIES
#  THEPEG_FIO_LIBRARIES

set(THEPEG_DIR "/home/maren/GEANT4/ThePEG-2.0.1-install")

find_path(THEPEG_INCLUDE_DIR ThePEG/Repository/EventGenerator.h
          HINTS $ENV{THEPEG_DIR}/include ${THEPEG_DIR}/include)

find_library(THEPEG_LIBRARIES NAMES ThePEG
             HINTS $ENV{THEPEG_DIR}/lib/ThePEG ${THEPEG_DIR}/lib/ThePEG)

get_filename_component(THEPEG_LIBRARY_DIR ${THEPEG_LIBRARIES} PATH)
set(THEPEG_INCLUDE_DIRS ${THEPEG_INCLUDE_DIR})

# handle the QUIETLY and REQUIRED arguments and set THEPEG_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(ThePEG DEFAULT_MSG THEPEG_INCLUDE_DIR THEPEG_LIBRARIES)

mark_as_advanced(THEPEG_FOUND THEPEG_INCLUDE_DIR THEPEG_LIBRARIES)
