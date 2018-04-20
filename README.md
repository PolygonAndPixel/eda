# Estimation of Distribution Algorithms
This package uses different kinds of EDA methods for usage in multimodal, non-convex, non-smooth likelihood (or "energy") surfaces. The interface of the minimizers is based on
[IceTray](http://software.icecube.wisc.edu/) used by [IceCube](https://icecube.wisc.edu/).
This package allows you to get a feeling how different algorithms work, how
many function evaluations they need and how efficient those are under different
circumstances and with different parameters. I focus on algorithms that need
only few function evaluations (and hence are very efficient), so they can be
used in tasks where the function evaluation is expensive (i.e. needs a lot of
computation time). A good estimation of the integral of the function is desirable
but not the main goal. Finding the parameters with the highest likelihood is
the main target here.

## Prerequisites
 - [OpenBlas](https://www.openblas.net/)
 - [Lapack](http://www.netlib.org/lapack/)
 - [Lapacke](http://www.netlib.org/lapack/lapacke.html)
 - [gfortran](https://gcc.gnu.org/fortran/) (or any other Fortran compiler. Just change the flag in the corresponding makefiles within the folder `src/Minimzer/*/`)
 - [boost](https://www.boost.org/) (any version that is not too old should do)
You may want to change the links in the Makefile at root for your libraries.

## Usage
Compile using `make release` at root. Go to `build/apps` and execute `eda`
with the first parameter the path to your minimizer configuration file
(examples can be found at `../../xml/Minimizer/`) and the second argument
is the path to your likelihood configuration file (examples can be found at
`../../xml/likelihood/`. You may delete all the configurations in that file
for the functions you don't want to use).

## Known Problems
PolyChord currently uses the unfortunate parameters for gaussian shells.
This leads to a very long run.
A similar problem occurs for MultiNest where it does not find the shells.
You may increase the number of live points in the configuration file to
overcome this problem.
