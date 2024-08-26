# PROJECT DESCRIPTION

dq_coeff is a package of FORTRAN source code to calculate the finite difference weighting coefficients for arbitrary sets of points.  The method uses a self-contained arbitrary precision library and Taylor series acceleration via hyper-dual numbers.  

## BUILDING

The code has been tested on linux systems using gfortan versions 7.3.0 and 10.2.1.

Edit the first four lines of the included Makefile to specify your fortran compiler, compiler flags and the location of the system lapack library (only needed for the flow around a cylinder and triangular annulus heat conduction examples for solving the coupled fluid dynamics/heat flow equations).

type "make" this will build the dq_coff routine, all its dependencies and seven example programs

## EXAMPLES

Seven example programs are included from the Advances in Engineering Software manuscript.  All the examples source code, input and output files are included in the examples subdirectory.  Note the six of the examples (all except the hexagonal finite difference solution for the triangular annulus heat flow example) use randomly generated sets of points, so if your compiler produces a different sequence of random numbers, your output will differ.  

The example subdirectories ending in \_vp denote the standard (variable precision) implementation, while subdirectories ending in \_native are for modified versions of the code that do not use variable precision numbers to represent the hyper-dual number components, and results from these versions of the code cannot guarantee accurate results.

**Creeping flow around a cylinder:**
cylinder_32bit_vp
./run_cylinder_32bit.zsh will output results to the cylinder_32bit_summary.txt file and the results are shown in Figure 4c of the Advances in Engineering Software manuscript

cylinder_64bit_vp
./run_cylinder_64bit.zsh will output results to the cylinder_64bit_summary.txt file and the results are shown in Figure 4b of the Advances in Engineering Software manuscript

**Poisson equation:**
poisson_64bit_vp
./run_64bit.zsh > poisson_64bit_summary.txt will output results to the poisson_64bit_summary.txt file and the results are shown in Figure 3b of the Advances in Engineering Software manuscript

poisson_64bit_native
./run_64bit.zsh > poisson_64bit_summary.txt will output results to the poisson_64bit_summary.txt file and the results are shown in Figure 3b of the Advances in Engineering Software manuscript

poisson_128bit_native
./run_128bit.zsh > poisson_128bit_summary.txt will output results to the poisson_128bit_summary.txt file and the results are shown in Figure 3b of the Advances in Engineering Software manuscript

**Heat condition in a triangular annulus:**
triangle_64bit_vp
./run_triangle.zsh will output results to the triangle_summary.txt file and the results are shown in Figure 5b of the Advances in Engineering Software manuscript

triangle_finite_difference
./run_fd.csh > fd_summary.txt will output results to the fd_summary.txt file and the results are shown in Figure B.1 of the Advances in Engineering Software manuscript

## USAGE

Examples of usage are provided in six of examples above (excluding the hexagonal finite difference solution for the triangular annulus heat flow example).  A complete description of the FORTRAN interface and pseudo-code examples is given in Appendix A of the Advances in Engineering Software manuscript.

## CITATION

Jason Roberts, Meshless weighting coefficients for arbitrary nodes: the efficient computation to machine precision using hyper-dual numbers, Advances in Engineering Software, 2024, doi:10.1016/j.advengsoft.2024.103753

