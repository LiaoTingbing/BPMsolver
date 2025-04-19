get_filename_component(VCPKG_IMPORT_PREFIX "${CMAKE_CURRENT_LIST_DIR}/../../" ABSOLUTE)
# - Config file for the Armadillo package
# It defines the following variables
#  ARMADILLO_INCLUDE_DIRS - include directories for Armadillo
#  ARMADILLO_LIBRARY_DIRS - library directories for Armadillo (normally not used!)
#  ARMADILLO_LIBRARIES    - libraries to link against

# Tell the user project where to find our headers and libraries
set(ARMADILLO_INCLUDE_DIRS "${VCPKG_IMPORT_PREFIX}/include")
set(ARMADILLO_LIBRARY_DIRS "${VCPKG_IMPORT_PREFIX}/lib")

# Our library dependencies (contains definitions for IMPORTED targets)
include(CMakeFindDependencyMacro)
find_dependency(LAPACK)
include("${CMAKE_CURRENT_LIST_DIR}/ArmadilloLibraryDepends.cmake")

# These are IMPORTED targets created by ArmadilloLibraryDepends.cmake
set(ARMADILLO_LIBRARIES armadillo)

