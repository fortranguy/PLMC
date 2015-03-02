cmake_minimum_required(VERSION 3.0)

project(
    PLMC
    Fortran CXX
)

include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/SelectBuildType.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/CMakeMacroForceAddFlags.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/SetCompilerFlags.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/SetArchiverAndLinkingTool.cmake)

set(executable_prefix plmc)
set(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../lib ${CMAKE_CURRENT_BINARY_DIR}/lib)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/lib)

find_package(
    jsonfortran-${CMAKE_Fortran_COMPILER_ID}
    1.0.0
    REQUIRED
)
include_directories(${jsonfortran_INCLUDE_DIRS})

file(
    GLOB_RECURSE
    test_source_files
    src/test_*.f90
)

foreach(test_source ${test_source_files})
    get_filename_component(test_name ${test_source} NAME_WE)
    string(REGEX REPLACE ^test_ ""  test_name ${test_name})
    add_executable(
        ${executable_prefix}_${test_name}
        ${test_source}
    )
    set_property(
        TARGET ${executable_prefix}_${test_name}
        PROPERTY LINKER_LANGUAGE Fortran
    )
    target_link_libraries(
        ${executable_prefix}_${test_name}
        jsonfortran-static
        plmc_library
    )
    INSTALL(TARGETS ${executable_prefix}_${test_name} DESTINATION plmc)
endforeach(test_source)

