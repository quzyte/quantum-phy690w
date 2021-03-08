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
 *    
 *    This program solves the time-independent Gross–Pitaevskii nonlinear 
 *    partial differential equation in one space dimension in a trap using the
 *    imaginary-time propagation. The Gross–Pitaevskii equation describes the
 *    properties of a dilute trapped Bose–Einstein condensate. The equation is
 *    solved using the split-step Crank–Nicolson method by discretizing space
 *    and time. The discretized equation is then propagated in imaginary time
 *    over small time steps. When convergence is achieved, the method has
 *    yielded the stationary solution of the problem.
 *    
 *    Description of variables used in the code:
 *    
 *    opt     - decides which rescaling of GP equation will be used
 *    par     - parameter for rescaling of GP equation
 *    psi     - array with the wave function values
 *    pot     - array with the values of the potential
 *    G0      - final nonlinearity
 *    norm    - wave function norm
 *    rms     - root mean square radius
 *    mu      - chemical potential
 *    en      - energy
 *    Nx      - number of discretization points in the x-direction
 *    x       - array with the space mesh values in the x-direction
 *    dx      - spatial discretization step in the x-direction
 *    dt      - time discretization step
 *    Npas    - number of subsequent iterations with the fixed nonlinearity G0
 *    Nrun    - number of final iterations with the fixed nonlinearity G0
 *    output  - output file with the summary of final values of all physical 
 *              quantities
 *    initout - output file with the initial wave function
 *    Npasout - output file with the wave function obtained after the 
 *              subsequent Npas iterations, with the fixed nonlinearity G0
 *    Nrunout - output file with the final wave function obtained after the 
 *              final Nrun iterations
 *    outstpx - discretization step in the x-direction used to save wave 
 *              functions
 */

#include "imagtime1d.h"

int main(int argc, char **argv) {
   FILE *out;
   FILE *file;
   long cnti;
   double norm, rms, mu, en;
   double *psi;
   double *cbeta;
   double *dpsix;
   double *tmpxi, *tmpxj;
   
   if((argc != 3) || (strcmp(*(argv + 1), "-p") != 0)) {
      fprintf(stderr, "Usage: %s -p <parameterfile> \n", *argv);
      exit(EXIT_FAILURE);
   }

   if(! cfg_init(argv[2])) {
      fprintf(stderr, "Wrong input parameter file.\n");
      exit(EXIT_FAILURE);
   }
   
   readpar();
   
   x = alloc_double_vector(Nx);
   x2 = alloc_double_vector(Nx);

   pot = alloc_double_vector(Nx);
   psi = alloc_double_vector(Nx);
   
   dpsix = alloc_double_vector(Nx);
   
   calphax = alloc_double_vector(Nx - 1);
   cbeta =  alloc_double_vector(Nx - 1);
   cgammax = alloc_double_vector(Nx - 1);
   
   tmpxi = alloc_double_vector(Nx);
   tmpxj = alloc_double_vector(Nx);

   if(output != NULL) out = fopen(output, "w");
   else out = stdout;

   fprintf(out, "OPTION = %d\n", opt);
   fprintf(out, "NX = %12li\n", Nx);
   fprintf(out, "NPAS = %10li   NRUN = %10li\n", Npas, Nrun);
   fprintf(out, "DX = %8le\n", dx);
   fprintf(out, "DT = %8le\n", dt);
   fprintf(out, "G0 = %8le\n\n", G0);
   fprintf(out, "        %12s   %12s   %12s   %12s   %12s\n", "norm", "mu", "en", "rms", "psi(0)");
   fflush(out);
   
   G = 0.;
   
   init(psi);
   gencoef();
   calcnorm(&norm, psi, tmpxi);
   calcmuen(&mu, &en, psi, dpsix, tmpxi, tmpxj);
   calcrms(&rms, psi, tmpxi);
   fprintf(out, "INIT    %8le   %8le   %8le   %8le   %8le\n", norm, mu / par, en / par, rms, psi[Nx2]);
   fflush(out);
   if(initout != NULL) {
      file = fopen(initout, "w");
      for(cnti = 0; cnti < Nx; cnti += outstpx)
         fprintf(file, "%8le %8le\n", x[cnti], psi[cnti]);
      fclose(file);
   }
   
   G = par * G0;
   
   for(cnti = 0; cnti < Npas; cnti ++) {
      calcnu(psi);
      calclux(psi, cbeta);
      calcnorm(&norm, psi, tmpxi);
   }
   calcmuen(&mu, &en, psi, dpsix, tmpxi, tmpxj);
   calcrms(&rms, psi, tmpxi);
   fprintf(out, "NPAS    %8le   %8le   %8le   %8le   %8le\n", norm, mu / par, en / par, rms, psi[Nx2]);
   fflush(out);
   if(Npasout != NULL) {
      file = fopen(Npasout, "w");
      for(cnti = 0; cnti < Nx; cnti += outstpx)
         fprintf(file, "%8le %8le\n", x[cnti], psi[cnti]);
      fclose(file);
   }
   
   for(cnti = 0; cnti < Nrun; cnti ++) {
      calcnu(psi);
      calclux(psi, cbeta);
      calcnorm(&norm, psi, tmpxi);
   }
   calcmuen(&mu, &en, psi, dpsix, tmpxi, tmpxj);
   calcrms(&rms, psi, tmpxi);
   fprintf(out, "NRUN    %8le   %8le   %8le   %8le   %8le\n", norm, mu / par, en / par, rms, psi[Nx2]);
   fflush(out);
   if(Nrunout != NULL) {
      file = fopen(Nrunout, "w");
      for(cnti = 0; cnti < Nx; cnti += outstpx)
         fprintf(file, "%8le %8le\n", x[cnti], psi[cnti]);
      fclose(file);
   }

   if(output != NULL) fclose(out);

   free_double_vector(x);
   free_double_vector(x2);

   free_double_vector(pot);
   free_double_vector(psi);

   free_double_vector(dpsix);

   free_double_vector(calphax);
   free_double_vector(cbeta);
   free_double_vector(cgammax);

   free_double_vector(tmpxi);
   free_double_vector(tmpxj);

   return(EXIT_SUCCESS);
}

/**
 *    Reading input parameters from the configuration file.
 */
void readpar(void) {
   char *cfg_tmp;
   
   if((cfg_tmp = cfg_read("OPTION")) == NULL) {
      fprintf(stderr, "OPTION is not defined in the configuration file\n");
      exit(EXIT_FAILURE);
   }
   opt = atol(cfg_tmp);
   
   if((cfg_tmp = cfg_read("G0")) == NULL) {
      fprintf(stderr, "G0 is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   G0 = atof(cfg_tmp);

   if((cfg_tmp = cfg_read("NX")) == NULL) {
      fprintf(stderr, "NX is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   Nx = atol(cfg_tmp);

   if((cfg_tmp = cfg_read("DX")) == NULL) {
      fprintf(stderr, "DX is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   dx = atof(cfg_tmp);
   
   if((cfg_tmp = cfg_read("DT")) == NULL) {
      fprintf(stderr, "DT is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   dt = atof(cfg_tmp);
   
   if((cfg_tmp = cfg_read("NPAS")) == NULL) {
      fprintf(stderr, "NPAS is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   Npas = atol(cfg_tmp);
   
   if((cfg_tmp = cfg_read("NRUN")) == NULL) {
      fprintf(stderr, "NRUN is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   Nrun = atol(cfg_tmp);
   
   output = cfg_read("OUTPUT");
   initout = cfg_read("INITOUT");
   Npasout = cfg_read("NPASOUT");
   Nrunout = cfg_read("NRUNOUT");
   
   if((initout != NULL) || (Npasout != NULL) || (Nrunout != NULL)) {
      if((cfg_tmp = cfg_read("OUTSTPX")) == NULL) {
         fprintf(stderr, "OUTSTPX is not defined in the configuration file.\n");
         exit(EXIT_FAILURE);
      }
      outstpx = atol(cfg_tmp);
   }
   
   return;
}

/**
 *    Initialization of the space mesh, the potential, and the initial wave
 *    function.
 *    psi - array with the wave function values
 */
void init(double *psi) {
   long cnti;
   double pi, cpsi1, cpsi2;
   double tmp;
   
   if (opt == 2) par = 2.;
   else par = 1.;
   
   Nx2 = Nx / 2;
   dx2 = dx * dx;
   
   pi = 4. * atan(1.);
   cpsi1 = sqrt(sqrt(pi));
   cpsi2 = cpsi1 * sqrt(sqrt(2.));
   
   for(cnti = 0; cnti < Nx; cnti ++) {
      x[cnti] = (cnti - Nx2) * dx;
      x2[cnti] = x[cnti] * x[cnti];
      if(opt == 3) {
         pot[cnti] = 0.25 * x2[cnti];
         tmp = exp(- 0.25 * x2[cnti]);
         psi[cnti] = tmp / cpsi2;
      } else {
         pot[cnti] = x2[cnti];
         tmp = exp(- 0.5 * x2[cnti]);
         psi[cnti] = tmp / cpsi1;
      }
   }
   
   return;
}

/**
 *    Crank-Nicolson scheme coefficients generation.
 */
void gencoef(void) {
   long cnti;
   
   Ax0 = 1. + dt / dx2;
   
   Ax0r = 1. - dt / dx2;
   
   Ax = - 0.5 * dt / dx2;
   
   calphax[Nx - 2] = 0.;
   cgammax[Nx - 2] = - 1. / Ax0;
   for (cnti = Nx - 2; cnti > 0; cnti --) {
      calphax[cnti - 1] = Ax * cgammax[cnti];
      cgammax[cnti - 1] = - 1. / (Ax0 + Ax * calphax[cnti - 1]);
   }
   
   return;
}

/**
 *    Calculation of the wave function norm and normalization.
 *    norm - wave function norm
 *    psi  - array with the wave function values
 *    tmpx - temporary array
 */
void calcnorm(double *norm, double *psi, double *tmpx) {
   long cnti;
   
   for(cnti = 0; cnti < Nx; cnti ++) {
      tmpx[cnti] = psi[cnti] * psi[cnti];
   }
   
   *norm = sqrt(simpint(dx, tmpx, Nx));
   
   for(cnti = 0; cnti < Nx; cnti ++) {
      psi[cnti] /= *norm;
   }
   
   return;
}

/**
 *    Calculation of the chemical potential and energy.
 *    mu    - chemical potential
 *    en    - energy
 *    psi   - array with the wave function values
 *    dpsix - temporary array
 *    tmpxi - temporary array
 *    tmpxj - temporary array
 */
void calcmuen(double *mu, double *en, double *psi, double *dpsix, double *tmpxi, double *tmpxj) {
   long cnti;
   double psi2, psi2lin, dpsi2;
   
   diff(dx, psi, dpsix, Nx);
   
   for(cnti = 0; cnti < Nx; cnti ++) {
      psi2 = psi[cnti] * psi[cnti];
      psi2lin = psi2 * G;
      dpsi2 = dpsix[cnti] * dpsix[cnti];
      tmpxi[cnti] = (pot[cnti] + psi2lin) * psi2 + dpsi2;
      tmpxj[cnti] = (pot[cnti] + 0.5 * psi2lin) * psi2 + dpsi2;
   }

   *mu = simpint(dx, tmpxi, Nx);
   *en = simpint(dx, tmpxj, Nx);
   
   return;
}

/**
 *    Calculation of the root mean square radius.
 *    rms  - root mean square radius
 *    psi  - array with the wave function values
 *    tmpx - temporary array
 */
void calcrms(double *rms, double *psi, double *tmpx) {
   long cnti;
   double psi2;
   
   for(cnti = 0; cnti < Nx; cnti ++) {
      psi2 = psi[cnti] * psi[cnti];
      tmpx[cnti] = x2[cnti] * psi2;
   }
   
   *rms = sqrt(simpint(dx, tmpx, Nx));
   
   return;
}

/**
 *    Time propagation with respect to H1 (part of the Hamiltonian without spatial
 *    derivatives).
 *    psi - array with the wave function values
 */
void calcnu(double *psi) {
   long cnti;
   double psi2, psi2lin, tmp;
   
   for(cnti = 0; cnti < Nx; cnti ++) {
      psi2 = psi[cnti] * psi[cnti];
      psi2lin = psi2 * G;
      tmp = dt * (pot[cnti] + psi2lin);
      psi[cnti] *= exp(- tmp);
   }
   
   return;
}

/**
 *    Time propagation with respect to H2 (x-part of the Laplacian).
 *    psi   - array with the wave function values
 *    cbeta - Crank-Nicolson scheme coefficients
 */
void calclux(double *psi, double *cbeta) {
   long cnti;
   double c;
   
   cbeta[Nx - 2] = psi[Nx - 1];
   for (cnti = Nx - 2; cnti > 0; cnti --) {
      c = - Ax * psi[cnti + 1] + Ax0r * psi[cnti] - Ax * psi[cnti - 1];
      cbeta[cnti - 1] =  cgammax[cnti] * (Ax * cbeta[cnti] - c);
   }
   psi[0] = 0.;
   for (cnti = 0; cnti < Nx - 2; cnti ++) {
      psi[cnti + 1] = calphax[cnti] * psi[cnti] + cbeta[cnti];
   }
   psi[Nx - 1] = 0.;
   
   return;
}
