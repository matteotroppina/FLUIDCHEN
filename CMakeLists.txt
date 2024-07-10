# Definition of the minimum required cmake Version
cmake_minimum_required(VERSION 3.12)
# Definition of the Project
# Later you can access the project variable like ${CFDLAB_SOURCE_DIR}
project(CFDLAB VERSION 1.0)

# Define all configuration options
#option(option1 "compile using this option" ON)

# Definition of the C++ Standard 
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CUDA_COMPILER /opt/nvidia/hpc_sdk/Linux_x86_64/23.3/compilers/bin/nvc++)
set(CMAKE_CUDA_STANDARD 17)

# Creating the executable of our project and the required dependencies
# the executable is called fluidchen
file(GLOB files src/*.cpp)
add_executable(fluidchen ${files})

# You can find package like
# find_package(MPI)
# Require a package
find_package(MPI REQUIRED)
# Find a package with different components e.g. BOOST
# find_package(Boost COMPONENTS filesystem REQUIRED)

# VTK Library
find_package(VTK REQUIRED)
message (STATUS "VTK_VERSION: ${VTK_VERSION}")
include(${VTK_USE_FILE})

# if you use external libraries you have to link them like
target_link_libraries(fluidchen PRIVATE MPI::MPI_CXX)
target_link_libraries(fluidchen PRIVATE ${VTK_LIBRARIES})

# Add include directory
target_include_directories(fluidchen PUBLIC include)


unset (CMAKE_CXX_FLAGS CACHE)

set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

if("${CMAKE_BUILD_TYPE}" STREQUAL "GPU")
    message(STATUS "Using GPU configuration")
    set(CMAKE_CXX_COMPILER ${CMAKE_CUDA_COMPILER})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -acc=gpu -Minfo=accel -fast -g -cuda")
else()
    message(STATUS "Using CPU configuration")
endif()

# If you write tests, you can include your subdirectory (in this case tests) as done here
# Testing
