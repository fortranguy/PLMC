
if (NOT CMAKE_BUILD_TYPE)
    set(
        CMAKE_BUILD_TYPE "Debug" CACHE STRING
        "Options are: Debug, Release or Fast. Default: Debug."
        FORCE
    )
endif (NOT CMAKE_BUILD_TYPE)
set_property(
    CACHE CMAKE_BUILD_TYPE PROPERTY
    STRINGS "Debug" "Release" "Fast"
)

