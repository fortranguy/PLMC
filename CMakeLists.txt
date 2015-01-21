cmake_minimum_required(VERSION 3.0)

project(
    Monte-Carlo_Canonical
    Fortran CXX
)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/SelectBuildType.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/CMakeMacroForceAddFlags.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/SetCompilerFlags.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/SetArchiverAndLinkingTool.cmake)

set(executable_prefix mc)
set(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/mod ${CMAKE_CURRENT_BINARY_DIR}/mod)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/mod)

find_package(
    jsonfortran-${CMAKE_Fortran_COMPILER_ID}
    1.0.0
    REQUIRED
)
include_directories(${jsonfortran_INCLUDE_DIRS})

find_package(Doxygen)
if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
endif(DOXYGEN_FOUND)

file(
    GLOB_RECURSE
    programs_source_files
    src/program_*.f90
)

foreach(program_source ${programs_source_files})
    get_filename_component(program_name ${program_source} NAME_WE)
    string(REGEX REPLACE ^program_ ""  program_name ${program_name})
    add_executable(
        ${executable_prefix}_${program_name}
        ${program_source}
    )
    set_property(
        TARGET ${executable_prefix}_${program_name}
        PROPERTY LINKER_LANGUAGE Fortran
    )
    target_link_libraries(
        ${executable_prefix}_${program_name}
        jsonfortran-static
        mc_meta_module
    )
    INSTALL(TARGETS ${executable_prefix}_${program_name} DESTINATION monte-carlo)
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

set(CPACK_PACKAGE_NAME "Monte-Carlo")
set(CPACK_PACKAGE_VENDOR "PLMC")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Monte-Carlo Simulation")
set(CPACK_PACKAGE_VERSION "1.0.0")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "/usr/local")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "salomon.chung@u-pe.frrm")

include(CPack)

