cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(vHLLE_SMASH_particle_sampler)
include(CTest)

# Tell cmake where to find modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")


# Set the relevant generic compiler flags (optimisation + warnings)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Wold-style-cast -Werror=switch -std=c++17")
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

set(SAMPLER_LIBRARIES "")

find_package(SMASH)
if(${SMASH_FOUND})
  set(SAMPLER_LIBRARIES ${SAMPLER_LIBRARIES} ${SMASH_LIBRARIES})
  include_directories(${SMASH_INCLUDE_DIR})
  if(DEFINED ENV{SMASH_DIR})
    include_directories($ENV{SMASH_DIR}/3rdparty)
  else()
    message(FATAL_ERROR "Environment variable SMASH_DIR unset. Unable to include virtest.")
  endif()
endif(${SMASH_FOUND})

find_package(ROOT 5.34)
if(ROOT_FOUND)
  set(SAMPLER_LIBRARIES ${SAMPLER_LIBRARIES} ${ROOT_LIBRARIES})
  include_directories(SYSTEM "${ROOT_INCLUDE_DIR}")
  #target_link_libraries(sampler aws-cpp-sdk-s3 aws-cpp-sdk-core)
endif()

target_link_libraries(sampler ${SAMPLER_LIBRARIES})

message(STATUS "BUILD_TESTING=${BUILD_TESTING}")

if(BUILD_TESTING)
  # library for unit tests
  #add_library(sampler_static STATIC $<TARGET_OBJECTS:objlib>)
  #set_target_properties(sampler_static PROPERTIES OUTPUT_NAME sampler)

  # Add the executable for the gen test
  add_executable(gen tests/gen.cpp)

  # Link the test executable with Virtest and any other necessary libraries
  target_link_libraries(gen objlib ${SAMPLER_LIBRARIES})

  # Add the test to ctest
  add_test(NAME gen_test COMMAND gen)
  #add_subdirectory(tests)
endif()