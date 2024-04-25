# Toy-Cxx-LAPACK
Assorted C++-ified LAPACK routines for reference.

Currently just enough to obtain the eigensystem for a symmetric real matrix.

Requires at least -std=c++20.
 
Refer to the original FORTRAN routines for documentation.
Plaintext link: https://netlib.org/lapack/explore-html/

Several routines were changed to facilitate a standalone implementation.

This is intended for reference purposes only, and not meant for production code.
The only guarantee here is that everything is 0-based, and should be easier to follow if writing -something new- for a different value type or template code.

* Thse are complete translations from FORTRAN to C++ with index-space refactoring from Reference-LAPACK.

* There are _NO_ dependencies outside what little of the standard library was used.
* It follows that this does _NOT_ depend on actual LAPACK of any version.
 
* Should compile without warning using the default options. Just run CMake and build.
