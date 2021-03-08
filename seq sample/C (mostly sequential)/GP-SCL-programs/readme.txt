# 
# GP-SCL codes are developed by:
# 
# Dusan Vudragovic, Ivana Vidanovic, Antun Balaz
# (Institute of Physics Belgrade, Serbia, http://www.scl.rs/)
# 
# Paulsamy Muruganandam
# (Bharathidasan University, Tamil Nadu, India)
# 
# Sadhan K. Adhikari
# (Sao Paulo State University, Brazil)
#    
# Public use and modification of this code are allowed provided that the
# following papers are cited:
# P. Muruganandam et al., Comput. Phys. Commun. 180 (2009) 1888;
# D. Vudragovic et al., Comput. Phys. Commun. 183 (2012) 2021.
# The authors would be grateful for all information and/or comments regarding
# the use of the code.
#

GP-SCL set of C codes solves the time-(in)dependent Gross-Pitaevskii nonlinear
partial differential equation in one, two, and three space dimensions in a
trap using imaginary-time and real-time propagation. The Gross-Pitaevskii
equation describes the properties of a dilute trapped Bose-Einstein
condensate. The equation is solved using the split-step Crank-Nicolson method
by discretizing space and time, as described in [1, 2]. The discretized
equation is then propagated in imaginary or real time over small time steps.

GP-SCL set of codes represents improved and threaded C version of previously
published Fortran AEDU codes, described in:
[1] P. Muruganandam, S. K. Adhikari, Comput. Phys. Commun. 180 (2009) 1888.
[2] D. Vudragovic, I. Vidanovic, A. Balaz, P. Muruganandam, S. K. Adhikari,
    Comput. Phys. Commun. 183 (2012) 2021.

DESCRIPTION OF GP-SCL CODE DISTRIBUTION

I) Source codes

Source codes are written in the C programming language and are located in the
src folder. This folder has the following structure:

 - src/imagtime1d code solves the imaginary-time Gross-Pitaevskii equation in
   one space dimension with an anisotropic trap.
   
 - src/imagtime2d code solves the imaginary-time Gross-Pitaevskii equation in
   two space dimensions with an anisotropic trap.
   
 - src/imagtime2d-th code is threaded (OpenMP parallelized) version of
   src/imagtime2d code.
   
 - src/imagtime3d code solves the imaginary-time Gross-Pitaevskii equation in
   three space dimensions with an anisotropic trap.
   
 - src/imagtime3d-th code is threaded (OpenMP parallelized) version of
   src/imagtime3d code.
   
 - src/imagtimeaxial code solves the imaginary-time Gross-Pitaevskii equation
   in three space dimensions with an axially-symmetric trap. Effectively, this
   is a two-dimensional problem.
   
 - src/imagtimeaxial-th code is threaded (OpenMP parallelized) version of
   src/imagtimeaxial code.
   
 - src/imagtimecir code solves the imaginary-time Gross-Pitaevskii equation in
   two space dimensions with a circularly-symmetric trap. Effectively this is
   one-dimensional problem.
   
 - src/imagtimesph code solves the imaginary-time Gross-Pitaevskii equation in
   three space dimensions with a spherically-symmetric trap. Effectively this
   is one-dimensional problem.
   
 - src/realtime1d code solves the real-time Gross-Pitaevskii equation in one
   space dimension with an anisotropic trap.
   
 - src/realtime2d code solves the real-time Gross-Pitaevskii equation in two
   space dimensions with an anisotropic trap.
   
 - src/realtime2d-th code is threaded (OpenMP parallelized) version of
   src/realtime2d code.
   
 - src/realtime3d code solves the real-time Gross-Pitaevskii equation in three
   space dimensions with an anisotropic trap.
   
 - src/realtime3d-th code is threaded (OpenMP parallelized) version of
   src/realtime3d code.
   
 - src/realtimeaxial code solves the real-time Gross-Pitaevskii equation in
   three space dimensions with an axially-symmetric trap. Effectively this is
   two-dimensional problem.
   
 - src/realtimeaxial-th code is threaded (OpenMP parallelized) version of
   src/realtimeaxial code.
   
 - src/realtimecir code solves the real-time Gross-Pitaevskii equation in two
   space dimensions with a circularly-symmetric trap. Effectively this is
   one-dimensional problem.
 
 - src/realtimesph code solves the real-time Gross-Pitaevskii equation in
   three space dimensions with a spherically-symmetric trap. Effectively this
   is one-dimensional problem.
   
 - src/util provides utility functions for parsing of configuration files,
   integration and differentiation, as well as allocation/deallocation of
   memory.
   
II) Input files

For each programme, a specific parameter input file with characteristic set of
options and their detailed descriptions is provided in the input folder. This
folder also contains a general parameter input file input/input, which
contains all available options and can be used with any of the GP-SCL codes.
The following parameter input file templates are available in the input
folder:

 - input/imagtime1d-input parameter file template for use with the imagtime1d
   programme.
   
 - input/imagtime2d-input parameter file template for use with the imagtime2d
   and imagtime2d-th programmes.
   
 - input/imagtime3d-input parameter file template for use with the imagtime3d
   and imagtime3d-th programmes.
   
 - input/imagtimeaxial-input parameter file template for use with the
   imagtimeaxial and imagtimeaxial-th programmes.
   
 - input/imagtimecir-input parameter file template for use with the
   imagtimecir programme.
   
 - input/imagtimesph-input parameter file template for use with the
   imagtimesph programme.
   
 - input/realtime1d-input parameter file template for use with the realtime1d
   programme.
 
 - input/realtime2d-input parameter file template for use with the realtime2d
   and realtime2d-th programmes.
 
 - input/realtime3d-input parameter file template for use with the realtime3d
   and realtime3d-th programmes.
 
 - input/realtimeaxial-input parameter file template for use with the
   realtimeaxial and realtimeaxial-th programmes.
 
 - input/realtimecir-input parameter file template for use with the
   realtimecir programme.
 
 - input/realtimesph-input parameter file template for use with the
   realtimesph programme.
 
 - input/input parameter file template for use with all programmes available
   in GP-SCL distribution.

III) Examples

The examples folder contains examples of parameter input files and matching
outputs for all programmes available in the GP-SCL distribution. The folder
has the following structure:

 - examples/imagtime1d contains an example of imagtime1d programme parameter
   input file and matching output.
 
 - examples/imagtime2d contains an example of imagtime2d programme parameter
   input file and matching output.
 
 - examples/imagtime2d-th contains an example of imagtime2d-th programme
   parameter input file and matching output.
 
 - examples/imagtime3d contains an example of imagtime3d programme parameter
   input file and matching output.
 
 - examples/imagtime3d-th contains an example of imagtime3d-th programme
   parameter input file and matching output.
 
 - examples/imagtimeaxial contains an example of imagtimeaxial programme
   parameter input file and matching output.
 
 - examples/imagtimeaxial-th contains an example of imagtimeaxial-th programme
   parameter input file and matching output.
 
 - examples/imagtimecir contains an example of imagtimecir programme parameter
   input file and matching output.
 
 - examples/imagtimesph contains an example of imagtimesph programme parameter
   input file and matching output.
 
 - examples/realtime1d contains an example of realtime1d programme parameter
   input file and matching output.
 
 - examples/realtime2d contains an example of realtime2d programme parameter
   input file and matching output.
 
 - examples/realtime2d-th contains an example of realtime2d-th programme
   parameter input file and matching output.
 
 - examples/realtime3d contains an example of realtime3d programme parameter
   input file and matching output.
 
 - examples/realtime3d-th contains an example of realtime3d-th programme
   parameter input file and matching output.
 
 - examples/realtimeaxial contains an example of realtimeaxial programme
   parameter input file and matching output.
 
 - examples/realtimeaxial-th contains an example of realtimeaxial-th programme
   parameter input file and matching output.
 
 - examples/realtimecir contains an example of realtimecir programme parameter
   input file and matching output.
 
 - examples/realtimesph contains an example of realtimesph programme parameter
   input file and matching output.

IV) Installation

The provided makefile allows compilation of the GP-SCL codes. In addition to a
general makefile, we also provide five specific makefiles, tailored for use
with GNU C compiler (makefile.gcc), Intel C compiler (makefile.icc), IBM XL C
compiler (makefile.xlc), PGI C compiler (makefile.pgcc), and Sun (former
Oracle) C compiler (makefile.suncc).

Default optimization options set in makefiles are -O5 (gcc), -fast (icc), -O5
(xlc), -fast (pgcc), and -O3 -fast (suncc), and can be adjusted if needed in
the makefile variable CFLAGS. In the case of xlc compiler, CFLAGS variable
contains also the flag -q64, specifying the host architecture (for 32-bit
architecture, use -q32 flag).

The use of the general makefile:

   make <target> [compiler=<compiler>]

where possible targets are:
   
   all, serial, openmp, clean, help

as well as code-specific targets, which compile only a specified code:

   imagtime1d, imagtime2d, imagtime2d-th, imagtime3d, imagtime3d-th,
   imagtimeaxial, imagtimeaxial-th, imagtimecir, imagtimesph, realtime1d,
   realtime2d, realtime2d-th, realtime3d, realtime3d-th, realtimeaxial,
   realtimeaxial-th, realtimecir, realtimesph

Possible <compiler> values are:

   gcc, icc, xlc, pgcc, suncc

and should be chosen according to the available compilers on the compile host.
If compiler is not specified on the command line, the default gcc compiler is
used.

The use of compiler specific makefiles (makefile.gcc, makefile.icc,
makefile.xlc, makefile.pgcc, makefile.suncc) is the following:

   make -f <makefile> <target>

where <makefile> is one of:

   makefile.gcc, makefile.icc, makefile.xlc, makefile.pgcc, makefile.suncc

and possible targets are:

   all, serial, openmp, clean, help

as well as code-specific targets, which compile only a specified code:

   imagtime1d, imagtime2d, imagtime2d-th, imagtime3d, imagtime3d-th,
   imagtimeaxial, imagtimeaxial-th, imagtimecir, imagtimesph, realtime1d,
   realtime2d, realtime2d-th, realtime3d, realtime3d-th, realtimeaxial,
   realtimeaxial-th, realtimecir, realtimesph

Examples of compiling:

1) Compile all GP-SCL codes with the GNU C compiler:
   
   make all
   
   or
   
   make all compiler=gcc
   
   or
   
   make -f makefile.gcc all
   
2) Compile imagtimeaxial GP-SCL code with the Intel C compiler (if
   available):
   
   make imagtimeaxial compiler=icc
   
   or
   
   make -f makefile.icc imagtimeaxial
   
3) Compile threaded (OpenMP parallel) realtime3d-th GP-SCL code with the GNU C
   compiler:
   
   make realtime3d-th
   
   or
   
   make realtime3d-th compiler=gcc
   
   or
   
   make -f makefile.gcc realtime3d-th
   
4) Compile all threaded (OpenMP parallel) GP-SCL codes with the IBM XL C
   compiler (if available):
   
   make openmp compiler=xlc
   
   or
   
   make -f makefile.xlc openmp
   
5) Compile all serial GP-SCL codes with the Oracle C compiler (if available):
   
   make serial compiler=suncc
   
   or
   
   make -f makefile.suncc serial
   
V) Running the compiled codes

To run any of the GP-SCL codes compiled with the make command, use:
   
   ./<codename> -p <parameterfile>
   
where <codename> is a name of the compiled executable, and <parameterfile> is
a parameter input file prepared by the user. Examples of parameter input
files are described in section II above, and are provided in the folder input,
and in the folder examples (together with matching outputs).

Examples:

1) Run imagtime2d compiled code with the parameter input file
input/imagtime2d-input:
   
   ./imagtime2d -p input/imagtime2d-input
   
2) Run threaded (OpenMP parallelized) realtime3d-th compiled code with
the parameter input file realtime3d-input in the local directory using 8 CPU
cores:
   
   export OMP_NUM_THREADS=8
   ./realtime3d-th -p realtime3d-input
   
3) Run threaded (OpenMP parallelized) realtimeaxial-th compiled code with the
self-prepared parameter input file axial-input in the local directory using 3
CPU cores:
   
   export OMP_NUM_THREADS=3
   ./realtimeaxial-th -p axial-input
