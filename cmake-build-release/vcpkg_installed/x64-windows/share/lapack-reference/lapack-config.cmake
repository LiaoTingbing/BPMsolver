# Compute locations from <prefix>/lib/cmake/lapack-<v>/<self>.cmake
get_filename_component(_LAPACK_SELF_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)

# Load lapack targets from the install tree if necessary.
set(_LAPACK_TARGET "lapack")
if(_LAPACK_TARGET AND NOT TARGET "${_LAPACK_TARGET}")
  include("${_LAPACK_SELF_DIR}/lapack-targets.cmake")
endif()
unset(_LAPACK_TARGET)

# Hint for project building against lapack
set(LAPACK_Fortran_COMPILER_ID "GNU")

# Report the blas and lapack raw or imported libraries.
if("ON")
  include(CMakeFindDependencyMacro)
  find_dependency(BLAS) # For current build type
  set(LAPACK_blas_LIBRARIES "${BLAS_LIBRARIES}")
else()
  set(LAPACK_blas_LIBRARIES "blas") # target carries link libraries
endif()
set(LAPACK_lapack_LIBRARIES "lapack")
set(LAPACK_LIBRARIES ${LAPACK_lapack_LIBRARIES}) # target carries link libraries

unset(_LAPACK_SELF_DIR)
