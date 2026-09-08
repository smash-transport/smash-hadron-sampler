add_library(objlib OBJECT)
add_executable(sampler main.cpp)
target_link_libraries(sampler PRIVATE objlib)
# Make the sampler executable available in the top-level build directory for easier execution
set_target_properties(sampler PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

find_package(SMASH REQUIRED CONFIG COMPONENTS Core)
target_link_libraries(sampler PRIVATE SMASH::smash)

find_package(ROOT 6.06 REQUIRED CONFIG COMPONENTS Hist Physics Tree)
target_link_libraries(sampler PRIVATE ROOT::Hist ROOT::Physics ROOT::Tree)

# Create object library for the linker which contains pre-compiled cpp files
target_sources(objlib PRIVATE gen.cpp mass_integration.cpp oscaroutput.cpp 
                              params.cpp tree.cpp spin.cpp vorticity.cpp)
set_property(TARGET objlib PROPERTY POSITION_INDEPENDENT_CODE 1)
target_link_libraries(objlib PUBLIC SMASH::smash ROOT::Hist ROOT::Physics ROOT::Tree)

target_include_directories(objlib PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                                         $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
                                         $<INSTALL_INTERFACE:include>)
