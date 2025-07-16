[TOC]

# Installation

The Kariba library can be installed using make.

## Prerequisites

The library has just one prerequisites: the GNU Scientific Library, GSL. This can often be installed with help of a package manager.


## Download

First, download or clone the repository.

You can download a zip file from https://github.com/antonpannekoek/kariba (click the green '<> Code' button, then select "Download ZIP"), or
clone the repository: on the command line, in a suitable directory:

```
  git clone https://github.com/antonpannekoek/kariba.git

  cd kariba
```

## Set up the installation

Edit `make.config` and set the variables as necessary. You can comment-out or leave blank the variables you don't need.

- `GSL`: the main path to your GNU Scientific Library (GSL)
  directory. If your system package manager installed the GSL, you can
  leave this blank. You can probably find this path with
  `pkgconf gsl --variable=prefix` or `pkg-config gsl --variable=prefix`.

- `CXX`: the C++ compiler. In particular if you use OpenMP, make sure the compiler supports this.

- `CXXFLAGS`: some compiler flags you may want to set (or unset) when building the library.

- `OPENMP`: the compiler and linker option when using OpenMP. Leave empty or commented-out when not building with OpenMP.



## Build the code

From the root directory, just run `make`.

This will call the Makefiles in the separate subdirectories: `src/`,
`examples/`, `examples/model/` and `tests/`. You can build
subdirectories separately by either going into the directory and
(re)running `make`, or from the root, run any of

- `make lib`

- `make examples`

- `make model`

- `make tests`

The default is to build everything.

You can speed up building the code with parallel builds, using `make -j`.
If you prefer to keep your system responsive during the build,
supply the `-j` option with the number of cores you want to use for
building, e.g. `make -j7` would leave one core for the system during
the build on an 8-core machine.

### Products

The build process produces a `src/libkariba.a` file that is the actual
Kariba library. This is the actual file that you would need to include
at the linking stage when building your own model.

In `src/kariba` are the header files. If you include a header file
like `#include <kariba/Particles.hpp>`, then provide the `src/`
directory for the include path to the compiler, with the `-I` option.


Other files are the three example files, `examples/corona`,
`examples/particles` and `examples/singlezone`.

In the `examples/model/` subdirectory, there is the `bhwrap` executable that
can serve as an example application.



## Compiling and running the test and examples

If you build everything from the root directory, the unit tests and
the examples should already have been built. Otherwise, enter the
respective directories, tests/ and examples/, and simply run `make`.

Once build, inside the respective directories, you can run all the unit tests with

```
./test_main
```

Hopefully all tests pass.

Run `./test_main -h` to see options to select specific unit tests.

----

For the examples, inside their directory, run them with

```
./corona
./particles
./singlezone
```

The singlezone example produces some output to the terminal:

```
Physical quantities:
Average electron Lorenz factor: 8038.47
Equipartition Ue/UB: 698.369
Electron power: 1.80399e+43
Magnetic power: 2.58315e+40
Proton power: 4.12068e+42
Total power: 2.21864e+43 erg s^-1, 2.73063e-05 Eddington
```

All the examples also write a bunch of output files in the output/
directory. You can check the output files for correctness by comparing
them with those in the output_comparison directory.


## Installing in a separate directory (optional)

As an optional step, you can install the library files, including the header files, into a separate directory. This could be `/usr/local/`, or something in your home directory, like `$HOME/sw` or `$HOME/.local`.

The install will install the static and shared library files into a `lib/` subdirectory, and the Kariba header files into a `include/kariba/` subdirectory.

Set the `PREFIX` variable in `make.config` to the base directory, then simply run

```
make install
```


## Documentation

If you like to build the documentation, run

```
make docs
```

from the root directory.

The build documentation is then available in `docs/html`. Open the `docs/html/index.html` file to start browsing the documentation.
