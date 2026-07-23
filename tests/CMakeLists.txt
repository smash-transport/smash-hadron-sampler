# Target for the unit-tests framework
add_library(vir-test INTERFACE)
target_include_directories(vir-test SYSTEM INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/virtest)

# Create a directory for test executables
set(TEST_OUTPUT_DIR ${CMAKE_BINARY_DIR}/tests)
file(MAKE_DIRECTORY ${TEST_OUTPUT_DIR})

# Set runtime output directory for test executables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${TEST_OUTPUT_DIR})

# Add the executable for the gen test
add_executable(gen gen.cpp)
target_link_libraries(gen PRIVATE objlib vir-test)
add_test(NAME gen_test COMMAND gen)

# Add the executable for the vorticity test
add_executable(vorticity vorticity.cpp)
target_link_libraries(vorticity PRIVATE objlib vir-test)
add_test(NAME vorticity_test COMMAND vorticity)

# Add the executable for the spin test
add_executable(spin spin.cpp)
target_link_libraries(spin PRIVATE objlib)
add_test(NAME spin_test COMMAND spin)
