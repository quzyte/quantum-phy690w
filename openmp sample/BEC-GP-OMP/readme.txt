#
# BEC-GP-OMP codes are developed and (c)opyright-ed by:
#
# Luis E. Young-S., Sadhan K. Adhikari
# (UNESP - Sao Paulo State University, Brazil)
#
# Paulsamy Muruganandam
# (Bharathidasan University, Tamil Nadu, India)
#
# Dusan Vudragovic, Antun Balaz
# (Scientific Computing Laboratory, Institute of Physics Belgrade, Serbia)
#
# Public use and modification of this code are allowed provided that the
# following three papers are cited:
#
# [1] L. E. Young-S. et al., Comput. Phys. Commun. 204 (2016) 209.
# [2] P. Muruganandam, S. K. Adhikari, Comput. Phys. Commun. 180 (2009) 1888.
# [3] D. Vudragovic et al., Comput. Phys. Commun. 183 (2012) 2021.
#
# The authors would be grateful for all information and/or comments
# regarding the use of the code.
#

BEC-GP-OMP set of Fortran and C codes solves the time-(in)dependent Gross-Pitaevskii
nonlinear partial differential equation for BECs with contact  interaction
in one, two, and three space dimensions in a trap using imaginary-time and real-time
propagation. The Gross-Pitaevskii equation describes the properties of a dilute trapped
Bose-Einstein condensate. The equation is solved using the split-step Crank-Nicolson method
by discretizing space and time, as described in [1]. The discretized equation is then
propagated in imaginary or real time over small time steps. Additional details are given in Refs. [1,2].

[1] P. Muruganandam, S. K. Adhikari, Comput. Phys. Commun. 180 (2009) 1888.
[2] D. Vudragovic, I. Vidanovic, A. Balaz, P. Muruganandam, S. K. Adhikari, Comput. Phys. Commun. 183 (2012) 2021.

DESCRIPTION OF BEC-GP-OMP CODE DISTRIBUTION

I) Source codes

Source codes are written in both the Fortran and C programming language and are located
in the src folder of the corresponding directory (f_programs, c_programs).

The Fortran src folder has the following structure:

 - src/imag1d.f90 code solves the imaginary-time Gross-Pitaevskii equation in
   one space dimension for the dynamics along x-axis in a harmonic trap.

 - src/imagcir.f90 code solves the imaginary-time Gross-Pitaevskii equation in
   two space dimensions for circular symmetry in a harmonic trap.

 - src/imagsph.f90 code solves the imaginary-time Gross-Pitaevskii equation in
   three space dimensions for spherical symmetry in a harmonic trap.

 - src/imag2d.f90 code solves the imaginary-time Gross-Pitaevskii equation in
   two space dimensions for the dynamics along x-y plane in an anisotropic harmonic trap.

 - src/imagaxi.f90 code solves the imaginary-time Gross-Pitaevskii equation in
   three space dimensions in an axially-symmetric harmonic trap.
   
 - src/imag3d.f90 code solves the imaginary-time Gross-Pitaevskii equation in
   three space dimensions in an anisotropic harmonic trap.

 - src/real1d.f90 code solves the real-time Gross-Pitaevskii equation in
   one space dimension for the dynamics along x-axis in a harmonic trap.

 - src/realcir.f90 code solves the real-time Gross-Pitaevskii equation in
   two space dimensions for circular symmetry in a harmonic trap.

 - src/realsph.f90 code solves the real-time Gross-Pitaevskii equation in
   three space dimensions for spherical symmetry in a harmonic trap.

 - src/real2d.f90 code solves the real-time Gross-Pitaevskii equation in
   two space dimensions for the dynamics along x-y plane in an anisotropic harmonic trap.

 - src/realaxi.f90 code solves the real-time Gross-Pitaevskii equation in
   three space dimensions in an axially-symmetric harmonic trap.
   
 - src/real3d.f90 code solves the real-time Gross-Pitaevskii equation in
   three space dimensions in an anisotropic harmonic trap.


The C src folder has the analogous structure:

 - src/imag1d/* code solves the imaginary-time Gross-Pitaevskii equation in
   one space dimension for the dynamics along x-axis in a harmonic trap.

 - src/imagcir/* code solves the imaginary-time Gross-Pitaevskii equation in
   two space dimensions for circular symmetry in a harmonic trap.

 - src/imagsph/* code solves the imaginary-time Gross-Pitaevskii equation in
   three space dimensions for spherical symmetry in a harmonic trap.

 - src/imag2d/* code solves the imaginary-time Gross-Pitaevskii equation in
   two space dimensions for the dynamics along x-y plane in an anisotropic harmonic trap.

 - src/imagaxi/* code solves the imaginary-time Gross-Pitaevskii equation in
   three space dimensions in an axially-symmetric harmonic trap.

 - src/imag3d/* code solves the imaginary-time Gross-Pitaevskii equation in
   three space dimensions in an anisotropic harmonic trap.

 - src/real1d/* code solves the real-time Gross-Pitaevskii equation in
   one space dimension for the dynamics along x-axis in a harmonic trap.

 - src/realcir/* code solves the real-time Gross-Pitaevskii equation in
   two space dimensions for circular symmetry in a harmonic trap.

 - src/realsph/* code solves the real-time Gross-Pitaevskii equation in
   three space dimensions for spherical symmetry in a harmonic trap.

 - src/real2d/* code solves the real-time Gross-Pitaevskii equation in
   two space dimensions for the dynamics along x-y plane in an anisotropic harmonic trap.

 - src/realaxi/* code solves the real-time Gross-Pitaevskii equation in
   three space dimensions in an axially-symmetric harmonic trap.

 - src/real3d/* code solves the real-time Gross-Pitaevskii equation in
   three space dimensions in an anisotropic harmonic trap.

 - srs/utils/*  provides utility functions for parsing of configuration files,
   integration and differentiation, as well as allocation/deallocation of
   memory.

II) Input parameters

Input parameters for Fortran codes are given within each program.

For each C program, a specific parameter input file has to be given on the command
line. Examples of input files with characteristic set of options and their detailed
descriptions is provided in the input folder.

III) Examples of matching outputs

The output folder in both f_programs and c_programs directory contains examples
of matching outputs for all programs and default inputs available in the BEC-GP-OMP
distribution. Some large density files are omitted to save space.

IV) Compilation

To compile the Fortran and C codes one can use the provided makefiles using:

   make <target> [compiler=<compiler>]

where possible targets are:
   
   all, clean, help

as well as code-specific targets, which compile only a specified code:

   imag1d, imagcir, imagsph, imag2d, imagaxi, imag3d, 
   real1d, realcir, realsph, real2d, realaxi, real3d, 

In the case of Fortran codes, the make command compiles all programs by default
with the Intel Fortran compiler. Possible <compiler> values are:

   ifort, gfort

where gfort is used for compiling with GFortran, and should be chosen according
to the available compilers on the compile host.

Examples compiling Fortran programs:

1) Compile all Fortran BEC-GP-OMP codes with the GFortran compiler:
   
   make all compiler=gfort
      
2) Compile imag3d Fortran code with the Intel Fortran compiler:
   
   make imag3d compiler=ifort 
   or   
   make imag3d

However, one can also use the instructions to compile Fortran programs that are
given within each source file.


For C programs, the provided makefile allows compilation of the BEC-GP-OMP codes using
either the GNU C compiler or Intel C compiler.

Default optimization options set in makefile are -O5 (gcc) and -fast (icc), and can be
adjusted if needed in the makefile variable CFLAGS.

If compiler is not specified on the command line, the default icc compiler is used.

For C codes possible <compiler> values are:

   icc, gcc

and should be chosen according to the available compilers on the compile host.

Examples compiling C programs:

1) Compile all C BEC-GP-OMP codes with the Intel C compiler:
   
   make all
   
   or
   
   make all compiler=icc

2) Compile all C BEC-GP-OMP codes with the GNU C compiler:

   make all compiler=gcc
      
3) Compile imag2d C code with the Intel C compiler:
   
   make imag2d compiler=icc
 
V) Running the compiled codes

When running the programs, by default, all available CPU cores will be used, i.e.,
the number of OpenMP threads used will be equal to the number of CPU cores available.
If you want to specify a different number of CPU cores to be used (e.g., in case
you want to run more than one program in paralell on the same computer), you should
define it by setting the shell variable OMP_NUM_THREADS to the desired number of
threads. For instance, if you want to use 4 threads, the following command should
be executed prior to runing the programs:

export OMP_NUM_THREADS=4     [if you are using bash shell]
setenv OMP_NUM_THREADS 4     [if you are using csh or tcsh shell]


V.A) Running Fortran programs:

Since the paramaters are specified within each program, to run the compiled
Fortran programs you just need to invoke the executable.

To run the Fortran programs use the executable:

   ./<codename>

for the imaginary time routines, where <codename> is a name of the compiled executable.

For the supplied realtime routines with NSTP = 0, use:

   ./<codename> < imag<n>d-den.txt

where <n> corresponds to the dimension of the problem.

Example of running a Fortran realtime program:

Run real3d compiled Fortran code using the input file to read imag3d-den.txt:
   
   ./real3d < imag3d-den.txt

Important note: if you encounter segmentation fault erros when running Fortran
programs on Linux, you should increase the stack size. Running the following
command should resolve the issue:

   ulimit -s unlimited          [if you are using bash shell]
   limit stacksize unlimited    [if you are using csh or tcsh shell]

V.B) Running C programs:

To run any of the C codes compiled with the make command, you need to use the
following syntax:
   
   ./<codename> -i <parameterfile>
   
where <codename> is a name of the compiled executable, and <parameterfile> is
a parameter input file prepared by the user. Examples of parameter input
files are provided in the folder input.

Matching output  of the principal output files are given in the folder output; very
large density output files are omitted to save space.

Example of running a C program:

Run imag3d compiled code with the parameter input file input/imag3d-input:
  
   ./imag3d -i input/imag3d-input

VI) Execution times of the supplied programs

The C programs are fully parallelized when compiled by Intel icc or GNU gcc compiler.
The Fortran programs are fully parallelized only when Intel ifort compiler is used.
They will compile with GNU gfortran, but will use only one CPU core when executed.
When Intel compiler is used, execution times of the Fortran and C programs will be
similar to those from Table 1 for 1 and 20 CPU cores (depending on the actual type
of processors used). The execution times of the C programs using the GNU gcc compiler
are close to those obtained by the Intel icc. The execution times of the Fortran programs
using the GNU gfortran compiler are generally larger than the times obtained by the
programs compiled by the Intel ifort compiler (and executed on a single CPU core).
Hence, the use of the Intel ifort compiler is recommended for the Fortran programs.
