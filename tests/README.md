# Kariba Library Unit Tests

This directory contains unit tests for the Kariba C++ library using the [doctest](https://github.com/doctest/doctest) testing framework.

## Test Files

- `test_constants.cpp` - Tests for physical constants and conversion factors
- `test_thermal.cpp` - Tests for the Thermal particle distribution class
- `test_powerlaw.cpp` - Tests for the Powerlaw particle distribution class  
- `test_bbody.cpp` - Tests for the BBody blackbody radiation class
- `test_cyclosyn.cpp` - Tests for the Cyclosyn synchrotron radiation class
- `test_main.cpp` - Main test runner
- `doctest.h` - The doctest framework header (single-file library)

## Building and Running Tests

### Using Makefile (Traditional)

From the tests directory:
```bash
cd tests
make                    # Build all tests
make run-tests         # Build and run all tests
make run-main          # Run main test suite only
make run-constants     # Run constants tests only
```

Individual commands:
```bash
make kariba_tests      # Build main test suite
make test_constants    # Build constants tests
./kariba_tests         # Run main tests
./test_constants       # Run constants tests
```

### Using CMake (Alternative)

From the project root directory:
```bash
mkdir build && cd build
cmake ..
cmake --build . --target kariba_tests
```

To run the tests:
```bash
# Run all tests
cmake --build . --target run_tests

# Or run executables directly
./tests/kariba_tests
./tests/test_constants
```

### Using CMake from tests directory

From the tests directory:
```bash
mkdir build && cd build
cmake ..
cmake --build .
./kariba_tests
./test_constants
```

## Test Coverage

The tests currently cover:

### Core Classes
- **Thermal**: Maxwell-Boltzmann particle distributions with temperature setting, normalization, and average quantity calculations
- **Powerlaw**: Non-thermal power-law particle distributions with spectral indices, cooling, and maximum momentum calculations
- **BBody**: Blackbody radiation with temperature/luminosity setting and Wien's law relationships
- **Cyclosyn**: Synchrotron radiation with magnetic field dependence and frequency calculations

### Utilities
- **Constants**: Physical constants, unit conversions, and derived quantities
- **Base class functionality**: Array management, particle counting, energy calculations

## Test Framework

The tests use doctest, which provides:
- `CHECK()` and `REQUIRE()` assertion macros
- `TEST_CASE()` and `SUBCASE()` organization
- Floating-point comparisons with `doctest::Approx()`
- Automatic test discovery and execution

## Adding New Tests

To add tests for additional Kariba classes:

1. Create a new test file (e.g., `test_newclass.cpp`)
2. Include the doctest header and the class to test
3. Write test cases following the existing patterns
4. Add the new file to `CMakeLists.txt` in the `TEST_SOURCES` list
5. Rebuild and run tests

## Test Philosophy

These tests focus on:
- **Basic functionality**: Constructor/destructor, parameter setting
- **Physical consistency**: Scaling relationships, unit conversions
- **Edge cases**: Extreme parameter values, boundary conditions
- **Numerical stability**: Reasonable outputs for typical inputs

The tests are designed to catch regressions and validate the core physics of the Kariba library components.