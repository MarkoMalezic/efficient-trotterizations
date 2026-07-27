cmake_minimum_required(VERSION 3.10)

project(MyProject)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(include)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# Ensure a default build type if none is specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

# Option to enable -O3 optimization
option(USE_O3 "Enable -O3 optimization instead of -O2 in Release mode" OFF)

# Source files for the main program
file(GLOB_RECURSE SRC_FILES src/*.cpp)

# Create the main executable
add_executable(main ${SRC_FILES})

# Link against libquadmath
target_link_libraries(main PRIVATE quadmath)

# Apply different compiler options based on the build type

# Release with default O2 optimizer or O3 if USE_O3=ON
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    if(USE_O3)
        message(STATUS "Using Release mode with -O3 optimization")
        target_compile_options(main PRIVATE -O3)
        target_link_options(main PRIVATE -O3)
    else()
        message(STATUS "Using Release mode with default -O2 optimization")
        target_compile_options(main PRIVATE -O2)
        target_link_options(main PRIVATE -O2)
    endif()
endif()

# Debug with Og optimizer, AddressSanitizer and gprof profiler
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Using Debug mode: Enabling AddressSanitizer, -Og, and gprof profiling")
    
    # Enable AddressSanitizer
    target_compile_options(main PRIVATE -fsanitize=address -fno-omit-frame-pointer -Og -g)
    target_link_options(main PRIVATE -fsanitize=address -fno-omit-frame-pointer -Og -g)

    # Enable gprof profiling
    target_compile_options(main PRIVATE -pg)
    target_link_options(main PRIVATE -pg)
endif()
