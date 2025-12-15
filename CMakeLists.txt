# Minimum cmake version this is tested on
cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

# The name, version and language of the project
project(
  SMASH-hadron-sampler
  VERSION "3.3"
  LANGUAGES CXX)
set(SAMPLER_VERSION "${SMASH-hadron-sampler_VERSION}")
# If the following set() command is not commented out, the codebase is marked as
# dirty (for development on the way to next release)
#set(SAMPLER_VERSION "${SAMPLER_VERSION}-next")

# Fail if cmake is called in the source directory
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR "You don't want to configure in the source directory!")
endif()

# Tell cmake where to find modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# Determine System and set the relevant generic compiler flags accordingly
# (optimisation + warnings)
set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Wold-style-cast -Werror=switch -std=c++17"
)
if(NOT
   (APPLE
    AND CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "^arm"
    AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang"
    AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "15"))
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})

# ------------------------------------------------------------------------------
# Retrieve all metadata for --version option
# ------------------------------------------------------------------------------
include(GetGitRevisionDescription)
find_package(Git)
if(GIT_FOUND)
  # message("git found: ${GIT_EXECUTABLE} in version ${GIT_VERSION_STRING}")
  set(GIT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/.git")
  if(EXISTS "${GIT_DIR}")
    execute_process(
      COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      OUTPUT_VARIABLE GIT_BRANCH)
    string(STRIP ${GIT_BRANCH} GIT_BRANCH)
  else()
    set(GIT_BRANCH "NOTFOUND")
  endif()
endif(GIT_FOUND)
# The version number is already set via the project() command, since a VERSION
# option is provided in it, and passed to the variable SAMPLER_VERSION via the
# set() command after it. The following checks if .git information about the
# branch and the commit hash is available.
git_describe(SAMPLER_VERSION_VERBOSE)
if(SAMPLER_VERSION_VERBOSE)
  set(SAMPLER_VERSION "${SAMPLER_VERSION_VERBOSE}")
else()
  set(SAMPLER_VERSION "${CMAKE_PROJECT_NAME}-${SAMPLER_VERSION}")
endif()

# Set a timestamp for the sampler compilation
string(TIMESTAMP BUILD_DATE)

include_directories("${CMAKE_CURRENT_BINARY_DIR}/src/include")
# configure a header file to pass some of the CMake settings to the source code
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/include/build_metadata.h.in"
               "${CMAKE_CURRENT_BINARY_DIR}/src/include/build_metadata.h")
# ------------------------------------------------------------------------------

# Collect all the cpp files
set(sampler_src src/gen.cpp src/oscaroutput.cpp src/params.cpp src/tree.cpp src/spin.cpp src/vorticity.cpp)

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
    message(
      FATAL_ERROR
        "Environment variable SMASH_DIR unset. Unable to include virtest.")
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

# enable standard CTest
include(CTest)

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
  #add_subdirectory(tests)
endif()
