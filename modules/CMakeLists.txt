cmake_minimum_required(VERSION 3.0)

project(
    PLMC_Library
    Fortran
)

include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/SelectBuildType.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/CMakeMacroForceAddFlags.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/SetCompilerFlags.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/SetArchiverAndLinkingTool.cmake)

find_package(
    jsonfortran-${CMAKE_Fortran_COMPILER_ID}
    4.2.0
    REQUIRED
)

file(
    GLOB_RECURSE
    lib_source_files
    micro/*.f90
    random/*.f90
    environment/*.f90
    mixture/*.f90
    short_interactions/*.f90
    long_interactions/*.f90
    changes/*.f90
    metropolis_algorithms/*.f90
    observables/*.f90
    writers/*.f90
    macro/*.f90
)

add_library(
    plmc_library
    STATIC
    ${lib_source_files}
)

include_directories(${jsonfortran_INCLUDE_DIRS})

set_property(
    TARGET plmc_library
    PROPERTY LINKER_LANGUAGE Fortran
)
target_link_libraries(
    plmc_library
    jsonfortran-static
)
