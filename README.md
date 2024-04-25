# Toy-Cxx-LAPACK
Assorted C++-ified LAPACK routines for reference.

Currently just enough to obtain the eigensystem for a symmetric real matrix.

Translations from FORTRAN to C++ with index-space refactoring from Reference-LAPACK.
Requires at least -std=c++20.
 
Refer to the original FORTRAN routines for documentation.
Plaintext link: https://netlib.org/lapack/explore-html/

Several routines were changed to facilitate a standalone implementation.

This is intended for reference purposes only, and not meant for production code.
The only guarantee here is that everything is 0-based, and should be easier to follow if building the same for a different value ty
