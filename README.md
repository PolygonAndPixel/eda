# Estimation of Distribution Algorithms
I am testing different kinds of EDA methods for usage in multimodal, non-convex, non-smooth likelihood (or "energy") surfaces. The interface of the minimizers is based on
[IceTray](http://software.icecube.wisc.edu/) used by [IceCube](https://icecube.wisc.edu/).
This package allows you to get a feeling how different algorithms work, how
many function evaluations they need and how efficient those are under different
circumstances and with different parameters. I focus on algorithms that need
only few function evaluations (and hence are very efficient), so they can be
used in tasks where the function evaluation is expensive (i.e. needs a lot of
computation time). A good estimation of the integral of the function is desirable
but not the main goal. Finding the parameters with the highest likelihood is
the main target here.

## Usage
Compile using `make release` at root. Go to `build/apps` and execute `eda` using
one of the following parameters:
 - poly (run [PolyChord](https://arxiv.org/abs/1506.00171) with every test function)
 - maps (run [MAPS](https://pdfs.semanticscholar.org/3261/1cd9eaa917d5bdcccd688799768afbd63579.pdf) with every test function)
 - sample (sample the whole space for every test function)

## Prerequisites
 - [OpenBlas](https://www.openblas.net/)
 - [Lapack](http://www.netlib.org/lapack/)
 - [Lapacke](http://www.netlib.org/lapack/lapacke.html)
 - [gfortran](https://gcc.gnu.org/fortran/) (or any other Fortran compiler. Just change the flag in the corresponding makefiles within the folder `src/Minimzer/*/`)
You may want to change the links in the Makefile at root for your libraries.

## Known Problems
PolyChord currently uses the same parameters for all functions. This leads to
a very long run on gaussian shells.
