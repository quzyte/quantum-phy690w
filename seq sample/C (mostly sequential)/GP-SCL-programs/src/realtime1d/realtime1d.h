/**   
 *    GP-SCL codes are developed by:
 *    
 *    Dusan Vudragovic, Ivana Vidanovic, Antun Balaz
 *    (Institute of Physics Belgrade, Serbia, http://www.scl.rs/)
 *    
 *    Paulsamy Muruganandam
 *    (Bharathidasan University, Tamil Nadu, India)
 *    
 *    Sadhan K. Adhikari
 *    (Sao Paulo State University, Brazil)
 *    
 *    Public use and modification of this code are allowed provided that the
 *    following papers are cited:
 *    P. Muruganandam et al., Comput. Phys. Commun. 180 (2009) 1888;
 *    D. Vudragovic et al., Comput. Phys. Commun. 183 (2012) 2021.
 *    The authors would be grateful for all information and/or comments
 *    regarding the use of the code.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

char *output, *initout, *rmsout, *Nstpout, *Npasout, *Nrunout;
long outstpx, outstpt;

int opt;
long Nx;
long Nx2;
long Nstp, Npas, Nrun;
double dx;
double dx2;
double dt;
double G0, Gpar, G;
double par;

double *x;
double *x2;
double *pot;

double complex Ax0, Ax0r, Ax;
double complex *calphax;
double complex *cgammax;

void readpar(void);
void init(double complex *);
void gencoef(void);
void calcnorm(double *, double complex *, double *);
void calcmuen(double *, double *, double complex *, double *, double *, double *);
void calcrms(double *, double complex *, double *);
void calcnu(double complex *);
void calclux(double complex *, double complex *);

extern double simpint(double, double *, long);
extern void diff(double, double *, double *, long);

extern double *alloc_double_vector(long);
extern void free_double_vector(double *);
extern double complex *alloc_complex_vector(long);
extern void free_complex_vector(double complex *);

extern int cfg_init(char *);
extern char *cfg_read(char *);
