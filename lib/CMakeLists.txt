cmake_minimum_required(VERSION 3.0)

project(
    PLMC_Library
    Fortran CXX
)

include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/SelectBuildType.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/CMakeMacroForceAddFlags.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/SetCompilerFlags.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/SetArchiverAndLinkingTool.cmake)

file(
    GLOB_RECURSE
    lib_source_files
    src/data_*.f90
    src/class_*.f90
)

add_library(
    plmc_library
    STATIC
    ${lib_source_files}
)

set_property(
    TARGET plmc_library
    PROPERTY LINKER_LANGUAGE Fortran
)
target_link_libraries(
    plmc_library
)
