# - Locate Herwig library
# in a directory defined via  HERWIG_ROOT_DIR or HERWIG_DIR environment variable
# Defines:
#
#  HERWIG_FOUND
#  HERWIG_INCLUDE_DIR
#  HERWIG_INCLUDE_DIRS (not cached)
#  HERWIG_LIBRARIES
#  HERWIG_FIO_LIBRARIES

set(HERWIG_DIR "/home/maren/GEANT4/Herwig-7.0.1-install")

find_path(HERWIG_INCLUDE_DIR Herwig/Decay/HwDecayHandler.h
          HINTS $ENV{HERWIG_DIR}/include ${HERWIG_DIR}/include)

find_library(HERWIG_LIBRARIES NAMES HwRunDirectories
             HINTS $ENV{HERWIG_DIR}/lib/Herwig ${HERWIG_DIR}/lib/Herwig)

get_filename_component(HERWIG_LIBRARY_DIR ${HERWIG_LIBRARIES} PATH)
set(HERWIG_INCLUDE_DIRS ${HERWIG_INCLUDE_DIR})

# handle the QUIETLY and REQUIRED arguments and set HERWIG_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Herwig DEFAULT_MSG HERWIG_INCLUDE_DIR HERWIG_LIBRARIES)

mark_as_advanced(HERWIG_FOUND HERWIG_INCLUDE_DIR HERWIG_LIBRARIES)
