# Definition of the minimum required cmake Version
cmake_minimum_required(VERSION 3.12)
# Definition of the Project
# Later you can access the project variable like ${CFDLAB_SOURCE_DIR}
project(CFDLAB VERSION 1.0)

# Define all configuration options
#option(option1 "compile using this option" ON)

## Set the C++ compiler to nvc++ GPU compiler
#set(CMAKE_CXX_COMPILER nvc++)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -acc")

# Definition of the C++ Standard 
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

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

# Add include directory
target_include_directories(fluidchen PUBLIC include)

set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

if(${CMAKE_BUILD_TYPE} STREQUAL "GPU")
    set(CMAKE_CXX_COMPILER nvc++)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -acc=gpu -Minfo=accel -fast -gpu=cc75")
endif()

# if you use external libraries you have to link them like
target_link_libraries(fluidchen PRIVATE MPI::MPI_CXX)
target_link_libraries(fluidchen PRIVATE ${VTK_LIBRARIES})

# If you write tests, you can include your subdirectory (in this case tests) as done here
# Testing
