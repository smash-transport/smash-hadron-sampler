# Minimum cmake version this is tested on
cmake_minimum_required(VERSION 3.28 FATAL_ERROR)

# The name, version and language of the project
project(
  SMASH-hadron-sampler
  VERSION "3.3"
  LANGUAGES CXX)
set(SAMPLER_VERSION "${SMASH-hadron-sampler_VERSION}")
# If the following set() command is not commented out, the codebase is marked as
# dirty (for development on the way to next release)
set(SAMPLER_VERSION "${SAMPLER_VERSION}-next")

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
find_package(Git)
if(GIT_FOUND)
  # message("git found: ${GIT_EXECUTABLE} in version ${GIT_VERSION_STRING}")
  set(GIT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/.git")
  if(EXISTS "${GIT_DIR}")
    execute_process(
      COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      OUTPUT_VARIABLE GIT_BRANCH)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} describe
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      OUTPUT_VARIABLE SAMPLER_VERSION_VERBOSE)
    string(STRIP ${GIT_BRANCH} GIT_BRANCH)
    string(STRIP ${SAMPLER_VERSION_VERBOSE} SAMPLER_VERSION_VERBOSE)
  else()
    set(GIT_BRANCH "NOTFOUND")
  endif()
else()
  set(GIT_BRANCH "NOTFOUND")
endif(GIT_FOUND)
# The version number is already set via the project() command, since a VERSION
# option is provided in it, and passed to the variable SAMPLER_VERSION via the
# set() command after it. The following checks if .git information about the
# branch and the commit hash is available.
if(SAMPLER_VERSION_VERBOSE)
  set(SAMPLER_VERSION "${SAMPLER_VERSION_VERBOSE}")
else()
  set(SAMPLER_VERSION "${CMAKE_PROJECT_NAME}-${SAMPLER_VERSION}")
endif()

# Set a timestamp for the sampler compilation
string(TIMESTAMP BUILD_DATE)

# configure a header file to pass some of the CMake settings to the source code
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/include/build_metadata.h.in"
               "${CMAKE_CURRENT_BINARY_DIR}/src/include/build_metadata.h")
# ------------------------------------------------------------------------------

# Add the subdirectory for the source code
add_subdirectory(src)

# enable standard CTest
include(CTest)
message(STATUS "BUILD_TESTING=${BUILD_TESTING}")
if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
