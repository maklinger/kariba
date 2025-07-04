cmake_minimum_required(VERSION 3.26)

# Set project name
project(KaribaTests)

# Find required packages
find_package(PkgConfig REQUIRED)
pkg_check_modules(GSL REQUIRED gsl)

find_package(OpenMP)

# Include configuration
if(EXISTS ${CMAKE_SOURCE_DIR}/../cmake.config)
    include(${CMAKE_SOURCE_DIR}/../cmake.config)
endif()

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../Kariba)

# Add compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall")

if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Collect Kariba source files
file(GLOB KARIBA_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/../Kariba/src/*.cpp")

# Test source files (excluding the constants test which has its own main)
set(TEST_SOURCES
    test_main.cpp
    test_thermal.cpp
    test_powerlaw.cpp
    test_bbody.cpp
    test_cyclosyn.cpp
)

# Create the main test executable
add_executable(kariba_tests ${TEST_SOURCES} ${KARIBA_SOURCES})

# Create separate executable for constants test (it has its own main)
add_executable(test_constants test_constants.cpp ${KARIBA_SOURCES})

# Link libraries
target_link_libraries(kariba_tests ${GSL_LIBRARIES})
target_link_libraries(test_constants ${GSL_LIBRARIES})

if(OpenMP_CXX_FOUND)
    target_link_libraries(kariba_tests OpenMP::OpenMP_CXX)
    target_link_libraries(test_constants OpenMP::OpenMP_CXX)
endif()

# Add library directories
target_link_directories(kariba_tests PRIVATE ${GSL_LIBRARY_DIRS})
target_link_directories(test_constants PRIVATE ${GSL_LIBRARY_DIRS})

# Add compile definitions
target_compile_definitions(kariba_tests PRIVATE ${GSL_CFLAGS_OTHER})
target_compile_definitions(test_constants PRIVATE ${GSL_CFLAGS_OTHER})

# Enable testing
enable_testing()

# Add tests
add_test(NAME KaribaConstantsTest COMMAND test_constants)
add_test(NAME KaribaMainTests COMMAND kariba_tests)

# Custom target to run all tests
add_custom_target(run_tests
    COMMAND ${CMAKE_CTEST_COMMAND} --verbose
    DEPENDS kariba_tests test_constants
    COMMENT "Running all Kariba unit tests"
)