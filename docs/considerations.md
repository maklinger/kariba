# Considerations

These are some considerations that were made during the restructuring
of BHJet into Kariba. It can also be seen as a form of a F.A.Q.

## Build system

The `make` build system is used. This is one of the simpler, and
readily available, build systems around (though it can get complex if
wanted). Given that the likely systems where Kariba should be build
and run on are Linux and macOS (from the terminal), no need for
another build system seems necessary.

Configuration is done with a small file where a few variables can be
set; there is no extensive configure (from autotools) script.

### CMake

By now, using CMake is far more standard for C++ projects. But CMake
can feel very verbose, with dozens of functions to set variables. It
is, in fact, a configuration language to create makefiles or
equivalent for a variety of OSes and build systems, which is not
necessary in this case: its output for Kariba would be a (more
complex) Makefile on the Linux and macOS shell.

Meson and Ninja are other alternatives for a build system, but the
availability and (in this case) the simplicity of using make directly
doesn't warrant the use of Meson or friends.

## Namespacing

Namespacing is liberally used, including prefixes to all common math
functions and iostream functions. This makes things a bit verbose, but
should also be clearer on what functions are actually used (even if
common).

A few exceptions still exist (including ones that slipped through the
replace-net), such as `size_t` instead of `std::size_t`.


### Constants

Constants have their own subnamespace. They are not capitalized, as is
standard for constants, since this is less mathematical, and also
prevents one from distinguishing prefixes such as `m` versus `M`
(although `MILLI` and `MEGA` could probably work fine here).


## Use of `size_t` and variable types

Many `int` occurrences have been replaced with `size_t`. This is not
out of necessity: array sizes are generally small enough that an `int`
will suffice. It is more indicative: the use of `size_t` indicates
this variable is used to limit the size of an array (although vectors
carry their own size in the form of a `.size()` method; which returns
a `size_T`) or is used to index an array (such as index variables in a
loop).

### Replacement of other variable types

Some other occurences of either `double` or `int` could be made more
explicit with fixed-width types, such as `float64_t`, `in16_t`, or
`int_least32_t`. This, however, seems unnecessary, as the C++ standard
guarantees enough precision and value-range for default doubles and
integers. Note that `float64_t` is only available since C++23.

### Use of `auto`

More use of `auto` is a possibility. This may make the code more
flexible. `auto` is probably the most useful with complex types, such
as `std::vector<double>`; but there are few places where the latter
can indeed be replaced by `auto`.

## Unit testing

[Doctest](https://github.com/doctest/doctest) is used for running unit
tests. The initial set of unit tests were created using Claude LLM
from the examples, and may have to be replaced by better
implementations over time.

Doctest is certainly not the definitive unit testing framework to be
used. Its advantage of being a single header file makes it easy to
include it with the repository (though its size and presumed
templating behind the scenes makes compilation relatively slow).


## Other improvements and issues

Check the issue tracker on the GitHub repository for suggestions of
improvements and features, and outstanding issues.
