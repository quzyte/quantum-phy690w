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
 *    partial differential equation in three space dimensions in a
 *    spherically-symmetric trap using the imaginary-time propagation. The
 *    Gross–Pitaevskii equation describes the properties of a dilute trapped
 *    Bose–Einstein condensate. The equation is solved using the split-step
 *    Crank–Nicolson method by discretizing space and time. The discretized
 *    equation is then propagated in imaginary time over small time steps.
 *    When convergence is achieved, the method has yielded the stationary
 *    solution of the problem.
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
 *    Nr      - number of discretization points in the r-direction
 *    r       - array with the space mesh values in the r-direction
 *    dr      - spatial discretization step in the r-direction
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
 *    outstpr - discretization step in the r-direction used to save wave 
 *              functions
 */

#include "imagtimesph.h"

int main(int argc, char **argv) {
   FILE *out;
   FILE *file;
   long cnti;
   double norm, rms, mu, en;
   double *psi;
   double *cbeta;
   double *dpsir;
   double *tmpri, *tmprj;
   
   if((argc != 3) || (strcmp(*(argv + 1), "-p") != 0)) {
      fprintf(stderr, "Usage: %s -p <parameterfile> \n", *argv);
      exit(EXIT_FAILURE);
   }

   if(! cfg_init(argv[2])) {
      fprintf(stderr, "Wrong input parameter file.\n");
      exit(EXIT_FAILURE);
   }
   
   readpar();
   
   r = alloc_double_vector(Nr);
   r2 = alloc_double_vector(Nr);

   pot = alloc_double_vector(Nr);
   psi = alloc_double_vector(Nr);
   
   dpsir = alloc_double_vector(Nr);
   
   calphar = alloc_double_vector(Nr - 1);
   cbeta =  alloc_double_vector(Nr - 1);
   cgammar = alloc_double_vector(Nr - 1);
   
   tmpri = alloc_double_vector(Nr);
   tmprj = alloc_double_vector(Nr);
   
   if(output != NULL) out = fopen(output, "w");
   else out = stdout;

   fprintf(out, "OPTION = %d\n", opt);
   fprintf(out, "NR = %12li\n", Nr);
   fprintf(out, "NPAS = %10li   NRUN = %10li\n", Npas, Nrun);
   fprintf(out, "DR = %8le\n", dr);
   fprintf(out, "DT = %8le\n", dt);
   fprintf(out, "G0 = %8le\n\n", G0);
   fprintf(out, "        %12s   %12s   %12s   %12s   %12s\n", "norm", "mu", "en", "rms", "psi(0)");
   fflush(out);
   
   G = 0.;
   
   init(psi);
   gencoef();
   calcnorm(&norm, psi, tmpri);
   calcmuen(&mu, &en, psi, dpsir, tmpri, tmprj);
   calcrms(&rms, psi, tmpri);
   fprintf(out, "INIT    %8le   %8le   %8le   %8le   %8le\n", norm, mu / par, en / par, rms, psi[1]/r[1]);
   fflush(out);
   if(initout != NULL) {
      file = fopen(initout, "w");
      fprintf(file, "%8le %8le\n", r[0], psi[1]/r[1]);
      for(cnti = 1; cnti < Nr; cnti += outstpr)
         fprintf(file, "%8le %8le\n", r[cnti], psi[cnti]/r[cnti]);
      fclose(file);
   }
   
   G = par * G0;
   
   for(cnti = 0; cnti < Npas; cnti ++) {
      calcnu(psi);
      calclur(psi, cbeta);
      calcnorm(&norm, psi, tmpri);
   }
   calcmuen(&mu, &en, psi, dpsir, tmpri, tmprj);
   calcrms(&rms, psi, tmpri);
   fprintf(out, "NPAS    %8le   %8le   %8le   %8le   %8le\n", norm, mu / par, en / par, rms, psi[1]/r[1]);
   fflush(out);
   if(Npasout != NULL) {
      file = fopen(Npasout, "w");
      fprintf(file, "%8le %8le\n", r[0], psi[1]/r[1]);
      for(cnti = 1; cnti < Nr; cnti += outstpr)
         fprintf(file, "%8le %8le\n", r[cnti], psi[cnti]/r[cnti]);
      fclose(file);
   }
   
   for(cnti = 0; cnti < Nrun; cnti ++) {
      calcnu(psi);
      calclur(psi, cbeta);
      calcnorm(&norm, psi, tmpri);
   }
   calcmuen(&mu, &en, psi, dpsir, tmpri, tmprj);
   calcrms(&rms, psi, tmpri);
   fprintf(out, "NRUN    %8le   %8le   %8le   %8le   %8le\n", norm, mu / par, en / par, rms, psi[1]/r[1]);
   fflush(out);
   if(Nrunout != NULL) {
      file = fopen(Nrunout, "w");
      fprintf(file, "%8le %8le\n", r[0], psi[1]/r[1]);
      for(cnti = 1; cnti < Nr; cnti += outstpr)
         fprintf(file, "%8le %8le\n", r[cnti], psi[cnti]/r[cnti]);
      fclose(file);
   }

   if(output != NULL) fclose(out);

   free_double_vector(r);
   free_double_vector(r2);

   free_double_vector(pot);
   free_double_vector(psi);

   free_double_vector(dpsir);

   free_double_vector(calphar);
   free_double_vector(cbeta);
   free_double_vector(cgammar);

   free_double_vector(tmpri);
   free_double_vector(tmprj);

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

   if((cfg_tmp = cfg_read("NR")) == NULL) {
      fprintf(stderr, "NR is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   Nr = atol(cfg_tmp);

   if((cfg_tmp = cfg_read("DR")) == NULL) {
      fprintf(stderr, "DR is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   dr = atof(cfg_tmp);
   
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
      if((cfg_tmp = cfg_read("OUTSTPR")) == NULL) {
         fprintf(stderr, "OUTSTPR is not defined in the configuration file.\n");
         exit(EXIT_FAILURE);
      }
      outstpr = atol(cfg_tmp);
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
   double cpsi1, cpsi2;
   double tmp;
   
   if (opt == 2) par = 2.;
   else par = 1.;
   
   Nr2 = Nr / 2;
   dr2 = dr * dr;
   
   pi = 4. * atan(1.);
   cpsi1 = sqrt(pi * sqrt(pi));
   cpsi2 = cpsi1 * sqrt(2. * sqrt(2.));
   
   for(cnti = 0; cnti < Nr; cnti ++) {
      r[cnti] = cnti * dr;
      r2[cnti] = r[cnti] * r[cnti];
      if(opt == 3) {
         pot[cnti] = 0.25 * r2[cnti];
         tmp = exp(- 0.25 * r2[cnti]);
         psi[cnti] = r[cnti] * tmp / cpsi2;
      } else {
         pot[cnti] = r2[cnti];
         tmp = exp(- 0.5 * r2[cnti]);
         psi[cnti] = r[cnti] * tmp / cpsi1;
      }
   }
   
   return;
}

/**
 *    Crank-Nicolson scheme coefficients generation.
 */
void gencoef(void) {
   long cnti;
   
   Ar0 = 1. + dt / dr2;
   
   Ar0r = 1. - dt / dr2;
   
   Ar = - 0.5 * dt / dr2;
   
   calphar[Nr - 2] = 0.;
   cgammar[Nr - 2] = - 1. / Ar0;
   for (cnti = Nr - 2; cnti > 0; cnti --) {
      calphar[cnti - 1] = Ar * cgammar[cnti];
      cgammar[cnti - 1] = - 1. / (Ar0 + Ar * calphar[cnti - 1]);
   }
   
   return;
}

/**
 *    Calculation of the wave function norm and normalization.
 *    norm - wave function norm
 *    psi  - array with the wave function values
 *    tmpr - temporary array
 */
void calcnorm(double *norm, double *psi, double *tmpr) {
   long cnti;
   
   for(cnti = 0; cnti < Nr; cnti ++) {
      tmpr[cnti] = 4 * pi * psi[cnti] * psi[cnti];
   }
   
   *norm = sqrt(simpint(dr, tmpr, Nr));
   
   for(cnti = 0; cnti < Nr; cnti ++) {
      psi[cnti] /= *norm;
   }
   
   return;
}

/**
 *    Calculation of the chemical potential and energy.
 *    mu    - chemical potential
 *    en    - energy
 *    psi   - array with the wave function values
 *    dpsir - temporary array
 *    tmpri - temporary array
 *    tmprj - temporary array
 */
void calcmuen(double *mu, double *en, double *psi, double *dpsir, double *tmpri, double *tmprj) {
   long cnti;
   double psi2, psi2lin, dpsi2;
   
   diff(dr, psi, dpsir, Nr);
   
   tmpri[0] = dpsir[1] * dpsir[1];
   tmprj[0] = tmpri[0];
   for(cnti = 1; cnti < Nr; cnti ++) {
      psi2 = psi[cnti] * psi[cnti];
      psi2lin = psi2 / r2[cnti] * G;
      dpsi2 = dpsir[cnti] * dpsir[cnti];
      tmpri[cnti] = (pot[cnti] + psi2lin) * psi2 + dpsi2;
      tmprj[cnti] = (pot[cnti] + 0.5 * psi2lin) * psi2 + dpsi2;
   }

   *mu = 4. * pi * simpint(dr, tmpri, Nr);
   *en = 4. * pi * simpint(dr, tmprj, Nr);
   
   return;
}

/**
 *    Calculation of the root mean square radius.
 *    rms  - root mean square radius
 *    psi  - array with the wave function values
 *    tmpr - temporary array
 */
void calcrms(double *rms, double *psi, double *tmpr) {
   long cnti;
   double psi2;
   
   for(cnti = 0; cnti < Nr; cnti ++) {
      psi2 = psi[cnti] * psi[cnti];
      tmpr[cnti] = 4. * pi *r2[cnti] * psi2;
   }
   
   *rms = sqrt(simpint(dr, tmpr, Nr));
   
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
   
   psi[0] = 0;
   for(cnti = 1; cnti < Nr; cnti ++) {
      psi2 = psi[cnti] * psi[cnti];
      psi2lin = G * psi2 / r2[cnti];
      tmp = dt * (pot[cnti] + psi2lin);
      psi[cnti] *= exp(- tmp);
   }
   
   return;
}

/**
 *    Time propagation with respect to H2 (r-part of the Laplacian).
 *    psi   - array with the wave function values
 *    cbeta - Crank-Nicolson scheme coefficients
 */
void calclur(double *psi, double *cbeta) {
   long cnti;
   double c;
   
   cbeta[Nr - 2] = psi[Nr - 1];
   for (cnti = Nr - 2; cnti > 0; cnti --) {
      c = - Ar * psi[cnti + 1] + Ar0r * psi[cnti] - Ar * psi[cnti - 1];
      cbeta[cnti - 1] =  cgammar[cnti] * (Ar * cbeta[cnti] - c);
   }
   psi[0] = 0.;
   for (cnti = 0; cnti < Nr - 1; cnti ++) {
      psi[cnti + 1] = calphar[cnti] * psi[cnti] + cbeta[cnti];
   }
   
   return;
}
