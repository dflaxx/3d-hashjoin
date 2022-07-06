# 3d-hashjoin
Source code for our paper:

Daniel Flachs, Magnus Müller, and Guido Moerkotte.
“The 3D Hash Join: Building On Non-Unique Join Attributes.”
12th Conference on Innovative Data Systems Research (CIDR 2022). January 9-12, 2022, Chaminade, USA.
https://www.cidrdb.org/cidr2022/papers/p18-flachs.pdf.

This repository contains two hash table implementations:
a regular chaining hash table and a nested/"3D" hash table as presented in the above paper.
Further, there is a small physical algebra with hash join implementations using the two hash table variants.
As the algebra is templated, most information (types for database tuples, hash functions etc.) must be available at compile time.
An example of how the algebra can be used can be found in `main_algebra_example.cc`.

`main_experiment1.cc` contains the key/foreign key experiment from the paper.
`main_experiment4.cc` contains the experiment with deferred unnesting from the paper.

The implementation uses C++20 constraints and concepts.


## Build

`make`

The build process can be configured by setting the following variables inside the `makefile` to either `0` or `1`:

- `DO_OPTIMIZE`: compiler optimization level (`-O0` vs. `-O3`)
- `DO_DEBUGINFO`: compiler debug info generation (`-g0` vs. `-g3`)
- `DO_ASSERTS`: enable/disable `assert(...)` statements
- `DO_PROFILE_GPROF`: enable/disable generation of profiling information for gprof
- `DO_SANITIZERS`: enable/disable address and undefined-behavior sanitizers


## Libraries

- [CLI11](https://github.com/CLIUtils/CLI11), v1.9.1: command line argument parsing.
  See lib/CLI11.hpp and lib/LICENSE for more info.
