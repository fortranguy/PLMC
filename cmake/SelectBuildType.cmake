
if (NOT CMAKE_BUILD_TYPE)
    set(
        build_type_message
        "Build type options are: Debug, Release or Fast. Default: Debug."
    )
    message(STATUS "${build_type_message}")
    set(
        CMAKE_BUILD_TYPE "Debug" CACHE STRING
        ${build_type_message}
        FORCE
    )
else (NOT CMAKE_BUILD_TYPE)
    message(STATUS "Build type is: ${CMAKE_BUILD_TYPE}")
endif (NOT CMAKE_BUILD_TYPE)
set_property(
    CACHE CMAKE_BUILD_TYPE PROPERTY
    STRINGS "Debug" "Release" "Fast"
)

