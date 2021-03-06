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
#include <omp.h>
#include <time.h>

#define MAX(a, b, c) (a > b) ? ((a > c) ? a : c) : ((b > c) ? b : c)
#define MAX_FILENAME_SIZE 256
#define RMS_ARRAY_SIZE     4

#define BOHR_RADIUS        5.2917720859e-11

char *output, *rmsout, *initout, *rmsout, *Nstpout, *Npasout, *Nrunout, *dynaout;
long outstpx, outstpy, outstpz, outstpt;

int opt;
long Na;
long Nx, Ny, Nz;
long Nx2, Ny2, Nz2;
long Nstp, Npas, Nrun;
double dx, dy, dz;
double dx2, dy2, dz2;
double dt;
double G0, Gpar, G;
double aho, as;
double vgamma, vnu, vlambda;
double par;
double pi;

double *x, *y, *z;
double *x2, *y2, *z2;
double ***pot;

double complex Ax0, Ay0, Az0, Ax0r, Ay0r, Az0r, Ax, Ay, Az;
double complex *calphax, *calphay, *calphaz;
double complex *cgammax, *cgammay, *cgammaz;

void readpar(void);
void init(double complex ***, double ***);
void gencoef(void);
void calcnorm(double *, double complex ***, double **, double **, double **);
void calcmuen(double *, double *, double complex ***, double ***, double ***, double ***, double **, double **, double **, double **, double **, double **);
void calcrms(double *, double complex ***, double **, double **, double **, double **, double **, double **, double **, double **, double **);
void calcnu(double complex ***);
void calclux(double complex ***, double complex **);
void calcluy(double complex ***, double complex **);
void calcluz(double complex ***, double complex **);
void outdenxyz(double complex ***, FILE *);
void outdenx(double complex ***, double *, double *, FILE *);
void outdeny(double complex ***, double *, double *, FILE *);
void outdenz(double complex ***, double *, double *, FILE *);
void outdenxy(double complex ***, double *, FILE *);
void outdenxz(double complex ***, double *, FILE *);
void outdenyz(double complex ***, double *, FILE *);
void outpsi2xz(double complex ***, FILE *);
void outpsi2xy(double complex ***, FILE *);

extern double simpint(double, double *, long);
extern void diff(double, double *, double *, long);

extern double *alloc_double_vector(long);
extern double complex *alloc_complex_vector(long);
extern double **alloc_double_matrix(long, long);
extern double complex **alloc_complex_matrix(long, long);
extern double ***alloc_double_tensor(long, long, long);
extern double complex ***alloc_complex_tensor(long, long, long);
extern void free_double_vector(double *);
extern void free_complex_vector(double complex *);
extern void free_double_matrix(double **);
extern void free_complex_matrix(double complex **);
extern void free_double_tensor(double ***);
extern void free_complex_tensor(double complex ***);

extern int cfg_init(char *);
extern char *cfg_read(char *);
