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
 *    This program solves the time-dependent Gross–Pitaevskii nonlinear
 *    partial differential equation in three space dimensions in an
 *    axially-symmetric trap using the real-time propagation. The
 *    Gross–Pitaevskii equation describes the properties of a dilute trapped
 *    Bose–Einstein condensate. The equation is solved using the split-step
 *    Crank–Nicolson method by discretizing space and time. The discretized
 *    equation is then propagated in real time over small time steps.
 *    
 *    Description of variables used in the code:
 *    
 *    opt       - decides which rescaling of GP equation will be used
 *    par       - parameter for rescaling of GP equation
 *    psi       - array with the wave function values
 *    pot       - array with the values of the potential
 *    G0        - final nonlinearity
 *    Gpar      - coefficient that multiplies nonlinear term in non-stationary
 *                problem during final Nrun iterations
 *    norm      - wave function norm
 *    rms       - root mean square radius
 *    mu        - chemical potential
 *    en        - energy
 *    Nrho      - number of discretization points in the rho-direction
 *    Nz        - number of discretization points in the z-direction
 *    rho       - array with the space mesh values in the rho-direction
 *    z         - array with the space mesh values in the z-direction
 *    drho      - spatial discretization step in the rho-direction
 *    dz        - spatial discretization step in the z-direction
 *    dt        - time discretization step
 *    kappa     - kappa coefficient of anisotropy of the trap
 *    lambda    - lambda coefficient of anisotropy of the trap
 *    Nstp      - number of subsequent iterations with fixed nonlinearity G0
 *    Npas      - number of subsequent iterations with the fixed nonlinearity
 *                G0
 *    Nrun      - number of final iterations with the fixed nonlinearity G0
 *    output    - output file with the summary of final values of all physical
 *                quantities
 *    initout   - output file with the initial wave function
 *    Nstpout   - output file with the wave function obtained after the first
 *                Nstp iterations
 *    Npasout   - output file with the wave function obtained after the
 *                subsequent Npas iterations, with the fixed nonlinearity G0
 *    Nrunout   - output file with the final wave function obtained after the
 *                final Nrun iterations
 *    rmsout    - output file with the time dependence of RMS during the final
 *                Nrun iterations
 *    outstprho - discretization step in the rho-direction used to save wave
 *                functions
 *    outstpz   - discretization step in the z-direction used to save wave
 *                functions
 *    outstpt   - time discretization step used to save RMS of the wave
 *                function
 */

#include "realtimeaxial-th.h"

int main(int argc, char **argv) {
   FILE *out;
   FILE *file;
   int nthreads;
   long cnti, cntj;
   double norm, rmsrho, rmsz, rms, mu, en;
   double complex **psi;
   double complex **cbeta;
   double **dpsirho, **dpsiz;
   double **tmprhoi, **tmpzi, **tmprhoj, **tmpzj;
   
   if((argc != 3) || (strcmp(*(argv + 1), "-p") != 0)) {
      fprintf(stderr, "Usage: %s -p <parameterfile> \n", *argv);
      exit(EXIT_FAILURE);
   }

   if(! cfg_init(argv[2])) {
      fprintf(stderr, "Wrong input parameter file.\n");
      exit(EXIT_FAILURE);
   }
   
   readpar();
   
   #pragma omp parallel
      #pragma omp master
         nthreads = omp_get_num_threads();
   
   rho = alloc_double_vector(Nrho);
   z = alloc_double_vector(Nz);

   rho2 = alloc_double_vector(Nrho);
   z2 = alloc_double_vector(Nz);

   pot = alloc_double_matrix(Nrho, Nz);
   psi = alloc_complex_matrix(Nrho, Nz);
   
   dpsirho = alloc_double_matrix(Nrho, Nz);
   dpsiz = alloc_double_matrix(Nrho, Nz);
   
   calpharho = alloc_complex_vector(Nrho - 1);
   calphaz = alloc_complex_vector(Nz - 1);
   cbeta =  alloc_complex_matrix(nthreads, MAX(Nrho, Nz) - 1);
   cgammarho = alloc_complex_vector(Nrho - 1);
   cgammaz = alloc_complex_vector(Nz - 1);
   
   tmprhoi = alloc_double_matrix(nthreads, Nrho);
   tmpzi = alloc_double_matrix(nthreads, Nz);
   tmprhoj = alloc_double_matrix(nthreads, Nrho);
   tmpzj = alloc_double_matrix(nthreads, Nz);
   
   if(output != NULL) out = fopen(output, "w");
   else out = stdout;
   
   fprintf(out, "OPTION = %d\n", opt);
   fprintf(out, "NRHO = %12li   NZ = %12li\n", Nrho, Nz);
   fprintf(out, "NSTP =   %10li   NPAS = %10li   NRUN = %10li\n", Nstp, Npas, Nrun);
   fprintf(out, "DRHO = %8le   DZ = %8le\n", drho, dz);
   fprintf(out, "DT = %8le\n", dt);
   fprintf(out, "AL = %8le   BL = %8le\n", kappa, lambda);
   fprintf(out, "G0 = %8le\n\n", G0);
   fprintf(out, "        %12s   %12s   %12s   %12s   %12s   %12s   %12s\n", "norm", "mu", "en", "rmsrho", "rmsz", "rms", "psi(0)");
   fflush(out);
   
   G = 0.;
   
   init(psi);
   gencoef();
   calcnorm(&norm, psi, tmprhoi, tmpzi);
   calcmuen(&mu, &en, psi, dpsirho, dpsiz, tmprhoi, tmpzi, tmprhoj, tmpzj);
   calcrms(&rmsrho, &rmsz, &rms, psi, tmprhoi, tmpzi, tmprhoj, tmpzj);
   fprintf(out, "INIT    %8le   %8le   %8le   %8le   %8le   %8le   %8le\n", norm, mu / par, en / par, rmsrho, rmsz, rms, cabs(psi[0][Nz2 + 1]));
   fflush(out);
   if(initout != NULL) {
      file = fopen(initout, "w");
      for(cnti = 0; cnti < Nrho; cnti += outstprho) {
         for(cntj = 0; cntj < Nz; cntj += outstpz) {
            fprintf(file, "%8le %8le %8le\n", rho[cnti], z[cntj], cabs(psi[cnti][cntj]));
         }
         fprintf(file, "\n");
      }
      fclose(file);
   }
   
   for(cnti = 0; cnti < Nstp; cnti ++) {
      G += par * G0 / Nstp;
      calcnu(psi);
      calclurho(psi, cbeta);
      calcluz(psi, cbeta);
   }
   calcnorm(&norm, psi, tmprhoi, tmpzi);
   calcmuen(&mu, &en, psi, dpsirho, dpsiz, tmprhoi, tmpzi, tmprhoj, tmpzj);
   calcrms(&rmsrho, &rmsz, &rms, psi, tmprhoi, tmpzi, tmprhoj, tmpzj);
   fprintf(out, "NSTP    %8le   %8le   %8le   %8le   %8le   %8le   %8le\n", norm, mu / par, en / par, rmsrho, rmsz, rms, cabs(psi[0][Nz2 + 1]));
   fflush(out);
   if(Nstpout != NULL) {
      file = fopen(Nstpout, "w");
      for(cnti = 0; cnti < Nrho; cnti += outstprho) {
         for(cntj = 0; cntj < Nz; cntj += outstpz) {
            fprintf(file, "%8le %8le %8le\n", rho[cnti], z[cntj], cabs(psi[cnti][cntj]));
         }
         fprintf(file, "\n");
      }
      fclose(file);
   }
   
   for(cnti = 0; cnti < Npas; cnti ++) {
      calcnu(psi);
      calclurho(psi, cbeta);
      calcluz(psi, cbeta);
   }
   calcnorm(&norm, psi, tmprhoi, tmpzi);
   calcmuen(&mu, &en, psi, dpsirho, dpsiz, tmprhoi, tmpzi, tmprhoj, tmpzj);
   calcrms(&rmsrho, &rmsz, &rms, psi, tmprhoi, tmpzi, tmprhoj, tmpzj);
   fprintf(out, "NPAS    %8le   %8le   %8le   %8le   %8le   %8le   %8le\n", norm, mu / par, en / par, rmsrho, rmsz, rms, cabs(psi[0][Nz2 + 1]));
   fflush(out);
   if(Npasout != NULL) {
      file = fopen(Npasout, "w");
      for(cnti = 0; cnti < Nrho; cnti += outstprho) {
         for(cntj = 0; cntj < Nz; cntj += outstpz) {
            fprintf(file, "%8le %8le %8le\n", rho[cnti], z[cntj], cabs(psi[cnti][cntj]));
         }
         fprintf(file, "\n");
      }
      fclose(file);
   }
   
   G *= Gpar;
   
   if(rmsout != NULL) file = fopen(rmsout, "w");
   for(cnti = 1; cnti <= Nrun; cnti ++) {
      calcnu(psi);
      calclurho(psi, cbeta);
      calcluz(psi, cbeta);
      if((rmsout != NULL) && (cnti % outstpt == 0)) {
         calcrms(&rmsrho, &rmsz, &rms, psi, tmprhoi, tmpzi, tmprhoj, tmpzj);
         fprintf(file, "%8le %8le %8le\n", cnti * dt * par, rmsrho, rmsz);
         fflush(file);
      }
   }
   if(rmsout != NULL) fclose(file);
   calcnorm(&norm, psi, tmprhoi, tmpzi);
   calcmuen(&mu, &en, psi, dpsirho, dpsiz, tmprhoi, tmpzi, tmprhoj, tmpzj);
   calcrms(&rmsrho, &rmsz, &rms, psi, tmprhoi, tmpzi, tmprhoj, tmpzj);
   fprintf(out, "NRUN    %8le   %8le   %8le   %8le   %8le   %8le   %8le\n", norm, mu / par, en / par, rmsrho, rmsz, rms, cabs(psi[0][Nz2 + 1]));
   fflush(out);
   if(Nrunout != NULL) {
      file = fopen(Nrunout, "w");
      for(cnti = 0; cnti < Nrho; cnti += outstprho) {
         for(cntj = 0; cntj < Nz; cntj += outstpz) {
            fprintf(file, "%8le %8le %8le\n", rho[cnti], z[cntj], cabs(psi[cnti][cntj]));
         }
         fprintf(file, "\n");
      }
      fclose(file);
   }
   
   if(output != NULL) fclose(out);
   
   free_double_vector(rho);
   free_double_vector(z);

   free_double_vector(rho2);
   free_double_vector(z2);

   free_double_matrix(pot);
   free_complex_matrix(psi);

   free_double_matrix(dpsirho);
   free_double_matrix(dpsiz);

   free_complex_vector(calpharho);
   free_complex_vector(calphaz);
   free_complex_matrix(cbeta);
   free_complex_vector(cgammarho);
   free_complex_vector(cgammaz);

   free_double_matrix(tmprhoi);
   free_double_matrix(tmpzi);
   free_double_matrix(tmprhoj);
   free_double_matrix(tmpzj);

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
   
   if((cfg_tmp = cfg_read("GPAR")) == NULL) {
      fprintf(stderr, "GPAR is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   Gpar = atof(cfg_tmp);

   if((cfg_tmp = cfg_read("NRHO")) == NULL) {
      fprintf(stderr, "NRHO is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   Nrho = atol(cfg_tmp);

   if((cfg_tmp = cfg_read("NZ")) == NULL) {
      fprintf(stderr, "NZ is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   Nz = atol(cfg_tmp);

   if((cfg_tmp = cfg_read("DRHO")) == NULL) {
      fprintf(stderr, "DRHO is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   drho = atof(cfg_tmp);

   if((cfg_tmp = cfg_read("DZ")) == NULL) {
      fprintf(stderr, "DZ is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   dz = atof(cfg_tmp);
   
   if((cfg_tmp = cfg_read("DT")) == NULL) {
      fprintf(stderr, "DT is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   dt = atof(cfg_tmp);
   
   if((cfg_tmp = cfg_read("AL")) == NULL) {
      fprintf(stderr, "AL is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   kappa = atof(cfg_tmp);
   
   if((cfg_tmp = cfg_read("BL")) == NULL) {
      fprintf(stderr, "BL is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   lambda = atof(cfg_tmp);

   if((cfg_tmp = cfg_read("NSTP")) == NULL) {
      fprintf(stderr, "NSTP is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   Nstp = atol(cfg_tmp);
   
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
   rmsout = cfg_read("RMSOUT");
   Nstpout = cfg_read("NSTPOUT");
   Npasout = cfg_read("NPASOUT");
   Nrunout = cfg_read("NRUNOUT");
   
   if((initout != NULL) || (Nstpout != NULL) || (Npasout != NULL) || (Nrunout != NULL)) {
      if((cfg_tmp = cfg_read("OUTSTPRHO")) == NULL) {
         fprintf(stderr, "OUTSTPRHO is not defined in the configuration file.\n");
         exit(EXIT_FAILURE);
      }
      outstprho = atol(cfg_tmp);
      
      if((cfg_tmp = cfg_read("OUTSTPZ")) == NULL) {
         fprintf(stderr, "OUTSTPZ is not defined in the configuration file.\n");
         exit(EXIT_FAILURE);
      }
      outstpz = atol(cfg_tmp);
   }
   
   if(rmsout != NULL) {
      if((cfg_tmp = cfg_read("OUTSTPT")) == NULL) {
         fprintf(stderr, "OUTSTPT is not defined in the configuration file.\n");
         exit(EXIT_FAILURE);
      }
      outstpt = atol(cfg_tmp);
   }
   
   return;
}

/**
 *    Initialization of the space mesh, the potential, and the initial wave
 *    function.
 *    psi - array with the wave function values
 */
void init(double complex **psi) {
   long cnti, cntj;
   double kappa2, lambda2;
   double pi3, cpsi1, cpsi2;
   double tmp;
   
   if (opt == 2) par = 2.;
   else par = 1.;
   
   kappa2 = kappa * kappa;
   lambda2 = lambda * lambda;
   
   Nrho2 = Nrho / 2; Nz2 = Nz / 2;
   drho2 = drho * drho; dz2 = dz * dz;
   
   pi = 4. * atan(1.);
   pi3 = pi * pi * pi;
   cpsi1 = sqrt(sqrt(pi3 / (lambda * kappa2)));
   cpsi2 = cpsi1 * sqrt(sqrt(8.));
   
   for(cnti = 0; cnti < Nrho; cnti ++) {
      rho[cnti] = cnti * drho;
      rho2[cnti] = rho[cnti] * rho[cnti];
   }
   
   for(cnti = 0; cnti < Nz; cnti ++) {
      z[cnti] = (cnti - Nz2) * dz;
      z2[cnti] = z[cnti] * z[cnti];
   }
   
   for(cnti = 0; cnti < Nrho; cnti ++) {
      for(cntj = 0; cntj < Nz; cntj ++) {
         if(opt == 3) {
            pot[cnti][cntj] = 0.25 * (kappa2 * rho2[cnti] + lambda2 * z2[cntj]);
            tmp = exp(- 0.25 * (kappa * rho2[cnti] + lambda * z2[cntj]));
            psi[cnti][cntj] = tmp / cpsi2;
         } else {
            pot[cnti][cntj] = (kappa2 * rho2[cnti] + lambda2 * z2[cntj]);
            tmp = exp(- 0.5 * (kappa * rho2[cnti] + lambda * z2[cntj]));
            psi[cnti][cntj] = tmp / cpsi1;
         }
      }
   }
   
   return;
}

/**
 *    Crank-Nicolson scheme coefficients generation.
 */
void gencoef(void) {
   long cnti;
   
   Arho0 = 1. + I * dt / drho2;
   Az0 = 1. + I * dt / dz2;
   
   Arho0r = 1. - I * dt / drho2;
   Az0r = 1. - I * dt / dz2;
   
   Arho = - 0.5 * I * dt / drho2;
   Az = - 0.5 * I * dt / dz2;
   
   dArho = 0.25 * I * dt / drho;
   
   calpharho[0] = 1.;
   cgammarho[0] = - 1. / (Arho0 + Arho + dArho / rho[1]);
   for(cnti = 1; cnti < Nrho - 1; cnti ++) {
      calpharho[cnti] = (Arho - dArho / rho[cnti]) * cgammarho[cnti - 1];
      cgammarho[cnti] = - 1. / (Arho0 + (Arho + dArho / rho[cnti + 1]) * calpharho[cnti]);
   }
   
   calphaz[Nz - 2] = 0.;
   cgammaz[Nz - 2] = - 1. / Az0;
   for (cnti = Nz - 2; cnti > 0; cnti --) {
      calphaz[cnti - 1] = Az * cgammaz[cnti];
      cgammaz[cnti - 1] = - 1. / (Az0 + Az * calphaz[cnti - 1]);
   }
   
   return;
}

/**
 *    Calculation of the wave function norm and normalization.
 *    norm   - wave function norm
 *    psi    - array with the wave function values
 *    tmprho - temporary array
 *    tmpz   - temporary array
 */
void calcnorm(double *norm, double complex **psi, double **tmprho, double **tmpz) {
   int threadid;
   long cnti, cntj;
   
   #pragma omp parallel private(threadid, cnti, cntj)
   {
      threadid = omp_get_thread_num();
      
      #pragma omp for
      for(cnti = 0; cnti < Nrho; cnti ++) {
         for(cntj = 0; cntj < Nz; cntj ++) {
            tmpz[threadid][cntj] = cabs(psi[cnti][cntj]);
            tmpz[threadid][cntj] *= tmpz[threadid][cntj];
         }
         tmprho[0][cnti] = rho[cnti] * simpint(dz, tmpz[threadid], Nz);
      }
      #pragma omp barrier
      
      #pragma omp single
      *norm = sqrt(2. * pi * simpint(drho, tmprho[0], Nrho));
      
      #pragma omp for
      for(cnti = 0; cnti < Nrho; cnti ++) {
         for(cntj = 0; cntj < Nz; cntj ++) {
            psi[cnti][cntj] /= *norm;
         }
      }
   }
   
   return;
}

/**
 *    Calculation of the chemical potential and energy.
 *    mu      - chemical potential
 *    en      - energy
 *    psi     - array with the wave function values
 *    dpsirho - temporary array
 *    dpsiz   - temporary array
 *    tmprhoi - temporary array
 *    tmpzi   - temporary array
 *    tmprhoj - temporary array
 *    tmpzj   - temporary array
 */
void calcmuen(double *mu, double *en, double complex **psi, double **dpsirho, double **dpsiz, double **tmprhoi, double **tmpzi, double **tmprhoj, double **tmpzj) {
   int threadid;
   long cnti, cntj;
   double psi2, psi2lin, dpsi2;
   
   #pragma omp parallel private(threadid, cnti, cntj, psi2, psi2lin, dpsi2)
   {
      threadid = omp_get_thread_num();
      
      #pragma omp for
      for(cntj = 0; cntj < Nz; cntj ++) {
         for(cnti = 0; cnti < Nrho; cnti ++) {
            tmprhoi[threadid][cnti] = cabs(psi[cnti][cntj]);
         }
         diff(drho, tmprhoi[threadid], tmprhoj[threadid], Nrho);
         for(cnti = 0; cnti < Nrho; cnti ++) {
            dpsirho[cnti][cntj] = tmprhoj[threadid][cnti];
         }
      }
      
      #pragma omp for
      for(cnti = 0; cnti < Nrho; cnti ++) {
         for(cntj = 0; cntj < Nz; cntj ++) {
            tmpzi[threadid][cntj] = cabs(psi[cnti][cntj]);
         }
         diff(dz, tmpzi[threadid], tmpzj[threadid], Nz);
         for(cntj = 0; cntj < Nz; cntj ++) {
            dpsiz[cnti][cntj] = tmpzj[threadid][cntj];
         }
      }
      #pragma omp barrier

      #pragma omp for
      for(cnti = 0; cnti < Nrho; cnti ++) {
         for(cntj = 0; cntj < Nz; cntj ++) {
            psi2 = cabs(psi[cnti][cntj]);
            psi2 *= psi2;
            psi2lin = psi2 * G;
            dpsi2 = dpsirho[cnti][cntj] * dpsirho[cnti][cntj] + 
                    dpsiz[cnti][cntj] * dpsiz[cnti][cntj];
            tmpzi[threadid][cntj] = (pot[cnti][cntj] + psi2lin) * psi2 + dpsi2;
            tmpzj[threadid][cntj] = (pot[cnti][cntj] + 0.5 * psi2lin) * psi2 + dpsi2;
         }
         tmprhoi[0][cnti] = rho[cnti] * simpint(dz, tmpzi[threadid], Nz);
         tmprhoj[0][cnti] = rho[cnti] * simpint(dz, tmpzj[threadid], Nz);
      }
   }
   
   *mu = 2. * pi * simpint(drho, tmprhoi[0], Nrho);
   *en = 2. * pi * simpint(drho, tmprhoj[0], Nrho);
   
   return;
}

/**
 *    Calculation of the root mean square radius.
 *    rms     - root mean square radius
 *    psi     - array with the wave function values
 *    tmprhoi - temporary array
 *    tmpzi   - temporary array
 *    tmprhoj - temporary array
 *    tmpzj   - temporary array
 */
void calcrms(double *rmsrho, double *rmsz, double *rms, double complex **psi, double **tmprhoi, double **tmpzi, double **tmprhoj, double **tmpzj) {
   int threadid;
   long cnti, cntj;
   double psi2;
   
   #pragma omp parallel private(threadid, cnti, cntj, psi2)
   {
      threadid = omp_get_thread_num();
      
      #pragma omp for
      for(cnti = 0; cnti < Nrho; cnti ++) {
         for(cntj = 0; cntj < Nz; cntj ++) {
            psi2 = cabs(psi[cnti][cntj]);
            psi2 *= psi2;
            tmpzi[threadid][cntj] = rho2[cnti] * psi2;
            tmpzj[threadid][cntj] = z2[cntj] * psi2;
         }
         tmprhoi[0][cnti] = rho[cnti] * simpint(dz, tmpzi[threadid], Nz);
         tmprhoj[0][cnti] = rho[cnti] * simpint(dz, tmpzj[threadid], Nz);
      }
   }
   
   *rmsrho = sqrt(2. * pi * simpint(drho, tmprhoi[0], Nrho));
   *rmsz = sqrt(2. * pi * simpint(drho, tmprhoj[0], Nrho));
   *rms = sqrt(*rmsrho * *rmsrho + *rmsz * *rmsz);
   
   return;
}

/**
 *    Time propagation with respect to H1 (part of the Hamiltonian without spatial
 *    derivatives).
 *    psi - array with the wave function values
 */
void calcnu(double complex **psi) {
   long cnti, cntj;
   double psi2, psi2lin, tmp;
   
   #pragma omp parallel for private(cnti, cntj, psi2, psi2lin, tmp)
   for(cnti = 0; cnti < Nrho; cnti ++) {
      for(cntj = 0; cntj < Nz; cntj ++) {
         psi2 = cabs(psi[cnti][cntj]);
         psi2 *= psi2;
         psi2lin = psi2 * G;
         tmp = dt * (pot[cnti][cntj] + psi2lin);
         psi[cnti][cntj] *= cexp(- I * tmp);
      }
   }
   
   return;
}

/**
 *    Time propagation with respect to H2 (rho-part of the Laplacian).
 *    psi   - array with the wave function values
 *    cbeta - Crank-Nicolson scheme coefficients
 */
void calclurho(double complex **psi, double complex **cbeta) {
   int threadid;
   long cnti, cntj;
   double complex c;
   
   #pragma omp parallel private(threadid, cnti, cntj, c)
   {
      threadid = omp_get_thread_num();
      
      #pragma omp for
      for(cntj = 0; cntj < Nz; cntj ++) {
         cbeta[threadid][0] = 0.;
         for(cnti = 1; cnti < Nrho - 1; cnti ++) {
            c = psi[cnti][cntj] - Arho * (psi[cnti + 1][cntj] - 2 * psi[cnti][cntj] + psi[cnti - 1][cntj]) + dArho / rho[cnti] * (psi[cnti + 1][cntj] - psi[cnti - 1][cntj]);
            cbeta[threadid][cnti] = cgammarho[cnti - 1] * ((Arho + dArho / rho[cnti]) * cbeta[threadid][cnti - 1] - c);
         }
         psi[Nrho - 1][cntj] = 0.;
         for (cnti = Nrho - 1; cnti > 0; cnti --) {
            psi[cnti - 1][cntj] = calpharho[cnti - 1] * psi[cnti][cntj] + cbeta[threadid][cnti - 1];
         }
      }
   }

   return;
}

/**
 *    Time propagation with respect to H3 (z-part of the Laplacian).
 *    psi   - array with the wave function values
 *    cbeta - Crank-Nicolson scheme coefficients
 */
void calcluz(double complex **psi, double complex **cbeta) {
   int threadid;
   long cnti, cntj;
   double complex c;
   
   #pragma omp parallel private(threadid, cnti, cntj, c)
   {
      threadid = omp_get_thread_num();
      
      #pragma omp for
      for(cnti = 0; cnti < Nrho; cnti ++) {
         cbeta[threadid][Nz - 2] = psi[cnti][Nz - 1];
         for (cntj = Nz - 2; cntj > 0; cntj --) {
            c = - Az * psi[cnti][cntj + 1] + Az0r * psi[cnti][cntj] - Az * psi[cnti][cntj - 1];
            cbeta[threadid][cntj - 1] =  cgammaz[cntj] * (Az * cbeta[threadid][cntj] - c);
         }
         for (cntj = 0; cntj < Nz - 2; cntj ++) {
            psi[cnti][cntj + 1] = calphaz[cntj] * psi[cnti][cntj] + cbeta[threadid][cntj];
         }
      }
   }
   
   return;
}
