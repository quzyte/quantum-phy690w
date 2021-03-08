/**   
 *    BEC-GP-OMP codes are developed and (c)opyright-ed by:
 *  
 *    Luis E. Young-S., Sadhan K. Adhikari
 *    (UNESP - Sao Paulo State University, Brazil)
 *  
 *    Paulsamy Muruganandam
 *    (Bharathidasan University, Tamil Nadu, India)
 *
 *    Dusan Vudragovic, Antun Balaz
 *    (Scientific Computing Laboratory, Institute of Physics Belgrade, Serbia)
 *
 *    Public use and modification of this code are allowed provided that the 
 *    following three papers are cited:
 *
 *    [1] L. E. Young-S. et al., Comput. Phys. Commun. 204 (2016) 209.
 *    [2] P. Muruganandam, S. K. Adhikari, Comput. Phys. Commun. 180 (2009) 1888.
 *    [3] D. Vudragovic et al., Comput. Phys. Commun. 183 (2012) 2021.
 *
 *    The authors would be grateful for all information and/or comments
 *    regarding the use of the code.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <omp.h>

#define MAX_FILENAME_SIZE     256
#define RMS_ARRAY_SIZE        1

#define BOHR_RADIUS           5.2917720859e-11

char *output, *dynaout, *initout, *Nstpout, *Npasout, *Nrunout;
long outstpx, outstpt;

int opt;
long Na;
long Nstp, Npas, Nrun;
long Nx;
long Nx2;
double g, g_3d, gpar, g_1d, drho, drho2;
double aho, as;
double dx;
double dx2;
double dt;
double vgamma;
double par;
double pi;

double *x;
double *x2;
double *pot;

double complex Ax0, Ax0r, Ax;
double complex *calphax;
double complex *cgammax;

void readpar(void);
void initpsi(double complex *psi, double *abc);
void initpot(void);
void gencoef(void);
void calcnorm(double *, double complex *, double *);
void calcmuen(double *, double *, double complex *, double *, double *, double *);
void calcrms(double *, double complex *, double *);
void calcnu(double complex *);
void calclux(double complex *, double complex *);
void outpsi2x(double complex *, FILE *);

extern double simpint(double, double *, long);
extern void diff(double, double *, double *, long);

extern int cfg_init(char *);
extern char *cfg_read(char *);

extern double *alloc_double_vector(long);
extern double complex *alloc_complex_vector(long);

extern void free_double_vector(double *);
extern void free_complex_vector(double complex *);

