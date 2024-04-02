PROJECT DESCRIPTION

dq_coeff is a package of FORTRAN source code to calculate the finite difference weighting coefficients for arbitrary sets of points.  The method uses a self-contained arbitrary precision library and Taylor series acceleration via hyper-dual numbers.  

BUILDING

The code has been tested on linux systems using gfortan versions 7.3.0 and 10.2.1.

Edit the first four lines of the included Makefile to specify your fortran compiler, compiler flags and the location of the system lapack library (only needed for the flow around a cylinder example for solving the coupled fluid dynamics equations).

type "make" this will build the dq_coff routine, all its dependencies and three example programs

EXAMPLES

Three example programs are included from the Journal of Scientific Computing manuscript.  All the examples source code, input and output files are included in the examples subdirectory.  Note the two of the examples (cylinder and poisson) use randomly generated sets of points, so if your compiler produces a different sequence of random numbers, your output will differ.  poisson_regular is provided as an example without random numbers.

Creeping flow around a cylinder: 

./cylinder < cylinder_input.txt 

Poisson equation: For a random set of points

./poisson  

While for a regular set of points

./poisson_regular

USAGE

Examples of usage are provided in the three examples above.  In general FORTRAN pseudo code is

use dq_coeff_mod ! routines to calculate weighting coefficients

use vector_mod   ! routines for accurate dot-product evaluation


! code to find N points of interest

call calc_dq_coeff(target_point,points,derivative,coeff,stat, &

    epsilon2_step,set_epsilon2,required_epsilon2,required_precision, &

    initial_epsilon2,initial_precision,taylor_series))


! INPUTS

! target_point : double(3) vector containing the (x,y,z) point of interest

! points : double (N,3) array containing the points being used to approximate the derivative at "target_point"

! derivative : vector of character strings containing M derivatives requested, can include "interp", "d_dx", "d_dy", "d_dz", "d2_dx2", "d2_dxdy", "d2_dxdz", "d2_dy2", "d2_dydz" and "d2_dz2"

! OUTPUTS

! coeff : double (N,M) array of finite difference weighting coefficients for M derivatives at N points

! stat : integer status, 0=no error, 1=required precision exceeded available arbitrary precision

! OPTIONAL INPUT

! epsilon2_step=double change from default factor of 4 when refining epsilon2

! set_epsilon2=double use a fixed value of epsilon2 without iteration to the asymptotic limit of epsilon2 -> infinity

! initial_epsilon2=double change from default starting value for epsilon2 iteration

! initial_precision=integer change from default starting arbitrary precision number of mantissa words for precision iteration

! taylor_series=integer choose fixed order of Taylor series expansion, value choices are 0 (no Taylor series), 1, 2 or 3.

! OPTIONAL OUTPUT

! required_epsilon2=double value of epsilon2 required for iterative convergence

! required_precision=integer number of mantissa words required for arbitrary precision iterative convergence


CITATION

Jason Roberts, Meshless weighting coefficients for arbitrary nodes: the efficient computation to machine precision using hyper-dual numbers, Advances in Engineering Software, (submitted)

