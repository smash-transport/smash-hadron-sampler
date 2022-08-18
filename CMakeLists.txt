cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(vHLLE_SMASH_particle_sampler)

# Tell cmake where to find modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")


# Set the relevant generic compiler flags (optimisation + warnings)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Wold-style-cast -Werror=switch -std=c++17 -mfpmath=sse")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})

# Collect all the cpp files
set(sampler_src
        src/gen.cpp
        src/params.cpp
        src/tree.cpp
        )

# Create object library for the linker which contains pre-compiled cpp files
add_library(objlib OBJECT ${sampler_src})
set_property(TARGET objlib PROPERTY POSITION_INDEPENDENT_CODE 1)

add_executable(sampler src/main.cpp $<TARGET_OBJECTS:objlib>)
include_directories(src/include)

find_package(SMASH)
if(${SMASH_FOUND})
  include_directories(${SMASH_INCLUDE_DIR})
  target_link_libraries(sampler ${SMASH_LIBRARIES})
endif(${SMASH_FOUND})

find_package(ROOT 5.34)
if(ROOT_FOUND)
  include_directories(SYSTEM "${ROOT_INCLUDE_DIR}")
  target_link_libraries(sampler ${ROOT_LIBRARIES} ${SMASH_LIBRARIES})
  #target_link_libraries(sampler aws-cpp-sdk-s3 aws-cpp-sdk-core)
endif()
