cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(vHLLE_SMASH_particle_sampler)
include(CTest)

# Tell cmake where to find modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")


# Determine System and set the relevant generic compiler flags accordingly (optimisation + warnings) 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Wold-style-cast -Werror=switch -std=c++17")
if(NOT (APPLE AND
   CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "^arm" AND
   CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" AND
   CMAKE_CXX_COMPILER_VERSION VERSION_LESS "15"))
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")  
endif()
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})  

# Collect all the cpp files
set(sampler_src
        src/gen.cpp
        src/params.cpp
        src/tree.cpp
        src/spin.cpp
        src/vorticity.cpp
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
else()
  message(FATAL_ERROR "SMASH could not be found.")
endif()


find_package(ROOT 5.34)
if(ROOT_FOUND)
  set(SAMPLER_LIBRARIES ${SAMPLER_LIBRARIES} ${ROOT_LIBRARIES})
  include_directories(SYSTEM "${ROOT_INCLUDE_DIR}")
  #target_link_libraries(sampler aws-cpp-sdk-s3 aws-cpp-sdk-core)
endif()

target_link_libraries(sampler ${SAMPLER_LIBRARIES})

message(STATUS "BUILD_TESTING=${BUILD_TESTING}")

if(BUILD_TESTING)
  # Create a directory for test executables
  set(TEST_OUTPUT_DIR ${CMAKE_BINARY_DIR}/tests)
  file(MAKE_DIRECTORY ${TEST_OUTPUT_DIR})

  # Set runtime output directory for test executables
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${TEST_OUTPUT_DIR})

  # Add the executable for the gen test
  add_executable(gen tests/gen.cpp)
  target_link_libraries(gen objlib ${SAMPLER_LIBRARIES})
  add_test(NAME gen_test COMMAND gen)

  # Add the executable for the vorticity test
  add_executable(vorticity tests/vorticity.cpp)
  target_link_libraries(vorticity objlib ${SAMPLER_LIBRARIES})
  add_test(NAME vorticity_test COMMAND vorticity)

  # Add the executable for the spin test
  add_executable(spin tests/spin.cpp)
  target_link_libraries(spin objlib ${SAMPLER_LIBRARIES})
  add_test(NAME spin_test COMMAND spin)
endif()
