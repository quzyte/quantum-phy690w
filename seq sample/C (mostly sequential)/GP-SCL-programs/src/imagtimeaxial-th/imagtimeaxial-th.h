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
#include <omp.h>

#define MAX(a, b) (a > b) ? a : b

char *output, *initout, *Npasout, *Nrunout;
long outstprho, outstpz;

int opt;
long Nrho, Nz;
long Nrho2, Nz2;
long Npas, Nrun;
double drho, dz;
double drho2, dz2;
double dt;
double G0, G;
double kappa, lambda;
double par;
double pi;

double *rho, *z;
double *rho2, *z2;
double **pot;

double Arho0, Az0, Arho0r, Az0r, Arho, Az, dArho;
double *calpharho, *calphaz;
double *cgammarho, *cgammaz;

void readpar(void);
void init(double **);
void gencoef(void);
void calcnorm(double *, double **, double **, double **);
void calcmuen(double *, double *, double **, double **, double **, double **, double **, double **, double **);
void calcrms(double *, double *, double *, double **, double **, double **, double **, double **);
void calcnu(double **);
void calclurho(double **, double **);
void calcluz(double **, double **);

extern double simpint(double, double *, long);
extern void diff(double, double *, double *, long);

extern double *alloc_double_vector(long);
extern double **alloc_double_matrix(long, long);
extern void free_double_vector(double *);
extern void free_double_matrix(double **);

extern int cfg_init(char *);
extern char *cfg_read(char *);
