cmake_minimum_required(VERSION 3.0)

project(PLMC Fortran)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/SelectBuildType.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/CMakeMacroForceAddFlags.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/SetCompilerFlags.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/SetArchiverAndLinkingTool.cmake)

set(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/modules ${CMAKE_CURRENT_BINARY_DIR}/lib)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/lib)

find_package(jsonfortran-${CMAKE_Fortran_COMPILER_ID} 5.0.0 REQUIRED)
include_directories(${jsonfortran_INCLUDE_DIRS})

file(GLOB_RECURSE programs_source_files programs/*.f90)
foreach(program_source ${programs_source_files})
    get_filename_component(program_name ${program_source} NAME_WE)
    add_executable(${program_name} ${program_source})
    set_property(TARGET ${program_name} PROPERTY LINKER_LANGUAGE Fortran)
    target_link_libraries(${program_name} jsonfortran-static plmc_library)
    INSTALL(TARGETS ${program_name} DESTINATION plmc)
endforeach(program_source)

# uninstall target
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/CMakeUninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/CMakeUninstall.cmake"
    IMMEDIATE @ONLY
)

add_custom_target(
    uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/CMakeUninstall.cmake
)

set(CPACK_PACKAGE_NAME "PLMC")
set(CPACK_PACKAGE_VENDOR "Laboratoire de Physique des Liquides et Milieux Complexes")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Physics of Liquids Monte-Carlo")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "1")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "/usr/local")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "salomon.chung@u-pe.fr")
include(CPack)
