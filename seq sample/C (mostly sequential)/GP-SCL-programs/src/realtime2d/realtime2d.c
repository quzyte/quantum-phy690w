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
 *    partial differential equation in two space dimensions in a trap using
 *    the real-time propagation. The Gross–Pitaevskii equation describes the
 *    properties of a dilute trapped Bose–Einstein condensate. The equation is
 *    solved using the split-step Crank–Nicolson method by discretizing space
 *    and time. The discretized equation is then propagated in real time over
 *    small time steps.
 *    
 *    Description of variables used in the code:
 *    
 *    opt     - decides which rescaling of GP equation will be used
 *    par     - parameter for rescaling of GP equation
 *    psi     - array with the wave function values
 *    pot     - array with the values of the potential
 *    G0      - final nonlinearity
 *    Gpar    - coefficient that multiplies nonlinear term in non-stationary
 *              problem during final Nrun iterations
 *    norm    - wave function norm
 *    rms     - root mean square radius
 *    mu      - chemical potential
 *    en      - energy
 *    Nx      - number of discretization points in the x-direction
 *    Ny      - number of discretization points in the y-direction
 *    x       - array with the space mesh values in the x-direction
 *    y       - array with the space mesh values in the y-direction
 *    dx      - spatial discretization step in the x-direction
 *    dy      - spatial discretization step in the y-direction
 *    dt      - time discretization step
 *    kappa   - kappa coefficient of anisotropy of the trap
 *    Nstp    - number of subsequent iterations with fixed nonlinearity G0
 *    Npas    - number of subsequent iterations with the fixed nonlinearity G0
 *    Nrun    - number of final iterations with the fixed nonlinearity G0
 *    output  - output file with the summary of final values of all physical
 *              quantities
 *    initout - output file with the initial wave function
 *    Nstpout - output file with the wave function obtained after the first
 *              Nstp iterations
 *    Npasout - output file with the wave function obtained after the
 *              subsequent Npas iterations, with the fixed nonlinearity G0
 *    Nrunout - output file with the final wave function obtained after the
 *              final Nrun iterations
 *    rmsout  - output file with the time dependence of RMS during the final
 *              Nrun iterations
 *    outstpx - discretization step in the x-direction used to save wave
 *              functions
 *    outstpy - discretization step in the y-direction used to save wave
 *              functions
 *    outstpt - time discretization step used to save RMS of the wave function
 */

#include "realtime2d.h"

int main(int argc, char **argv) {
   FILE *out;
   FILE *file;
   long cnti, cntj;
   double norm, rms, mu, en;
   double complex **psi;
   double complex *cbeta;
   double **dpsix, **dpsiy;
   double *tmpxi, *tmpyi, *tmpxj, *tmpyj;
   
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
   y = alloc_double_vector(Ny);

   x2 = alloc_double_vector(Nx);
   y2 = alloc_double_vector(Ny);

   pot = alloc_double_matrix(Nx, Ny);
   psi = alloc_complex_matrix(Nx, Ny);
   
   dpsix = alloc_double_matrix(Nx, Ny);
   dpsiy = alloc_double_matrix(Nx, Ny);
   
   calphax = alloc_complex_vector(Nx - 1);
   calphay = alloc_complex_vector(Ny - 1);
   cbeta =  alloc_complex_vector(MAX(Nx, Ny) - 1);
   cgammax = alloc_complex_vector(Nx - 1);
   cgammay = alloc_complex_vector(Ny - 1);
   
   tmpxi = alloc_double_vector(Nx);
   tmpyi = alloc_double_vector(Ny);
   tmpxj = alloc_double_vector(Nx);
   tmpyj = alloc_double_vector(Ny);
   
   if(output != NULL) out = fopen(output, "w");
   else out = stdout;

   fprintf(out, "OPTION = %d\n", opt);
   fprintf(out, "NX = %12li   NY = %12li\n", Nx, Ny);
   fprintf(out, "NSTP = %10li   NPAS = %10li   NRUN = %10li\n", Nstp, Npas, Nrun);
   fprintf(out, "DX = %8le   DY = %8le\n", dx, dy);
   fprintf(out, "DT = %8le\n", dt);
   fprintf(out, "AL = %8le\n", kappa);
   fprintf(out, "G0 = %8le\n\n", G0);
   fprintf(out, "        %12s   %12s   %12s   %12s   %12s\n", "norm", "mu", "en", "rms", "psi(0,0)");
   fflush(out);
   
   G = 0.;
   
   init(psi);
   gencoef();
   calcnorm(&norm, psi, tmpxi, tmpyi);
   calcmuen(&mu, &en, psi, dpsix, dpsiy, tmpxi, tmpyi, tmpxj, tmpyj);
   calcrms(&rms, psi, tmpxi, tmpyi);
   fprintf(out, "INIT    %8le   %8le   %8le   %8le   %8le\n", norm, mu / par, en / par, rms, cabs(psi[Nx2][Ny2]));
   fflush(out);
   if(initout != NULL) {
      file = fopen(initout, "w");
      for(cnti = 0; cnti < Nx; cnti += outstpx) {
         for(cntj = 0; cntj < Ny; cntj += outstpy) {
            fprintf(file, "%8le %8le %8le\n", x[cnti], y[cntj], cabs(psi[cnti][cntj]));
         }
         fprintf(file, "\n");
      }
      fclose(file);
   }
   
   for(cnti = 0; cnti < Nstp; cnti ++) {
      G += par * G0 / Nstp;
      calcnu(psi);
      calclux(psi, cbeta);
      calcluy(psi, cbeta);
   }
   calcnorm(&norm, psi, tmpxi, tmpyi);
   calcmuen(&mu, &en, psi, dpsix, dpsiy, tmpxi, tmpyi, tmpxj, tmpyj);
   calcrms(&rms, psi, tmpxi, tmpyi);
   fprintf(out, "NSTP    %8le   %8le   %8le   %8le   %8le\n", norm, mu / par, en / par, rms, cabs(psi[Nx2][Ny2]));
   fflush(out);
   if(Nstpout != NULL) {
      file = fopen(Nstpout, "w");
      for(cnti = 0; cnti < Nx; cnti += outstpx) {
         for(cntj = 0; cntj < Ny; cntj += outstpy) {
            fprintf(file, "%8le %8le %8le\n", x[cnti], y[cntj], cabs(psi[cnti][cntj]));
         }
         fprintf(file, "\n");
      }
      fclose(file);
   }
   
   for(cnti = 0; cnti < Npas; cnti ++) {
      calcnu(psi);
      calclux(psi, cbeta);
      calcluy(psi, cbeta);
   }
   calcnorm(&norm, psi, tmpxi, tmpyi);
   calcmuen(&mu, &en, psi, dpsix, dpsiy, tmpxi, tmpyi, tmpxj, tmpyj);
   calcrms(&rms, psi, tmpxi, tmpyi);
   fprintf(out, "NPAS    %8le   %8le   %8le   %8le   %8le\n", norm, mu / par, en / par, rms, cabs(psi[Nx2][Ny2]));
   fflush(out);
   if(Npasout != NULL) {
      file = fopen(Npasout, "w");
      for(cnti = 0; cnti < Nx; cnti += outstpx) {
         for(cntj = 0; cntj < Ny; cntj += outstpy) {
            fprintf(file, "%8le %8le %8le\n", x[cnti], y[cntj], cabs(psi[cnti][cntj]));
         }
         fprintf(file, "\n");
      }
      fclose(file);
   }
   
   G *= Gpar;
   
   if(rmsout != NULL) file = fopen(rmsout, "w");
   for(cnti = 1; cnti <= Nrun; cnti ++) {
      calcnu(psi);
      calclux(psi, cbeta);
      calcluy(psi, cbeta);
      if((rmsout != NULL) && (cnti % outstpt == 0)) {
         calcrms(&rms, psi, tmpxi, tmpyi);
         fprintf(file, "%8le %8le\n", cnti * dt * par, rms);
         fflush(file);
      }
   }
   if(rmsout != NULL) fclose(file);
   calcnorm(&norm, psi, tmpxi, tmpyi);
   calcmuen(&mu, &en, psi, dpsix, dpsiy, tmpxi, tmpyi, tmpxj, tmpyj);
   calcrms(&rms, psi, tmpxi, tmpyi);
   fprintf(out, "NRUN    %8le   %8le   %8le   %8le   %8le\n", norm, mu / par, en / par, rms, cabs(psi[Nx2][Ny2]));
   fflush(out);
   if(Nrunout != NULL) {
      file = fopen(Nrunout, "w");
      for(cnti = 0; cnti < Nx; cnti += outstpx) {
         for(cntj = 0; cntj < Ny; cntj += outstpy) {
            fprintf(file, "%8le %8le %8le\n", x[cnti], y[cntj], cabs(psi[cnti][cntj]));
         }
         fprintf(file, "\n");
      }
      fclose(file);
   }

   if(output != NULL) fclose(out);

   free_double_vector(x);
   free_double_vector(y);

   free_double_vector(x2);
   free_double_vector(y2);

   free_double_matrix(pot);
   free_complex_matrix(psi);

   free_double_matrix(dpsix);
   free_double_matrix(dpsiy);

   free_complex_vector(calphax);
   free_complex_vector(calphay);
   free_complex_vector(cbeta);
   free_complex_vector(cgammax);
   free_complex_vector(cgammay);

   free_double_vector(tmpxi);
   free_double_vector(tmpyi);
   free_double_vector(tmpxj);
   free_double_vector(tmpyj);

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

   if((cfg_tmp = cfg_read("NX")) == NULL) {
      fprintf(stderr, "NX is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   Nx = atol(cfg_tmp);

   if((cfg_tmp = cfg_read("NY")) == NULL) {
      fprintf(stderr, "NY is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   Ny = atol(cfg_tmp);

   if((cfg_tmp = cfg_read("DX")) == NULL) {
      fprintf(stderr, "DX is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   dx = atof(cfg_tmp);

   if((cfg_tmp = cfg_read("DY")) == NULL) {
      fprintf(stderr, "DY is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   dy = atof(cfg_tmp);
   
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
      if((cfg_tmp = cfg_read("OUTSTPX")) == NULL) {
         fprintf(stderr, "OUTSTPX is not defined in the configuration file.\n");
         exit(EXIT_FAILURE);
      }
      outstpx = atol(cfg_tmp);
      
      if((cfg_tmp = cfg_read("OUTSTPY")) == NULL) {
         fprintf(stderr, "OUTSTPY is not defined in the configuration file.\n");
         exit(EXIT_FAILURE);
      }
      outstpy = atol(cfg_tmp);
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
   double kappa2;
   double pi, cpsi1, cpsi2;
   double tmp;
   
   if (opt == 2) par = 2.;
   else par = 1.;
   
   kappa2 = kappa * kappa;
   
   Nx2 = Nx / 2; Ny2 = Ny / 2;
   dx2 = dx * dx; dy2 = dy * dy;
   
   pi = 4. * atan(1.);
   cpsi1 = sqrt(pi * sqrt(1. / kappa));
   cpsi2 = cpsi1 * sqrt(2.);
   
   for(cnti = 0; cnti < Nx; cnti ++) {
      x[cnti] = (cnti - Nx2) * dx;
      x2[cnti] = x[cnti] * x[cnti];
      for(cntj = 0; cntj < Ny; cntj ++) {
         y[cntj] = (cntj - Ny2) * dy;
         y2[cntj] = y[cntj] * y[cntj];
         if(opt == 3) {
            pot[cnti][cntj] = 0.25 * (x2[cnti] + kappa2 * y2[cntj]);
            tmp = exp(- 0.25 * (x2[cnti] + kappa * y2[cntj]));
            psi[cnti][cntj] = tmp / cpsi2;
         } else {
            pot[cnti][cntj] = x2[cnti] + kappa2 * y2[cntj];
            tmp = exp(- 0.5 * (x2[cnti] + kappa * y2[cntj]));
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
   
   Ax0 = 1. + I * dt / dx2;
   Ay0 = 1. + I * dt / dy2;
   
   Ax0r = 1. - I * dt / dx2;
   Ay0r = 1. - I * dt / dy2;
   
   Ax = - 0.5 * I * dt / dx2;
   Ay = - 0.5 * I * dt / dy2;
   
   calphax[Nx - 2] = 0.;
   cgammax[Nx - 2] = - 1. / Ax0;
   for (cnti = Nx - 2; cnti > 0; cnti --) {
      calphax[cnti - 1] = Ax * cgammax[cnti];
      cgammax[cnti - 1] = - 1. / (Ax0 + Ax * calphax[cnti - 1]);
   }
   
   calphay[Ny - 2] = 0.;
   cgammay[Ny - 2] = - 1. / Ay0;
   for (cnti = Ny - 2; cnti > 0; cnti --) {
      calphay[cnti - 1] = Ay * cgammay[cnti];
      cgammay[cnti - 1] = - 1. / (Ay0 + Ay * calphay[cnti - 1]);
   }
   
   return;
}

/**
 *    Calculation of the wave function norm and normalization.
 *    norm - wave function norm
 *    psi  - array with the wave function values
 *    tmpx - temporary array
 *    tmpy - temporary array
 */
void calcnorm(double *norm, double complex **psi, double *tmpx, double *tmpy) {
   long cnti, cntj;
   
   for(cnti = 0; cnti < Nx; cnti ++) {
      for(cntj = 0; cntj < Ny; cntj ++) {
         tmpy[cntj] = cabs(psi[cnti][cntj]);
         tmpy[cntj] *= tmpy[cntj];
      }
      tmpx[cnti] = simpint(dy, tmpy, Ny);
   }
   
   *norm = sqrt(simpint(dx, tmpx, Nx));
   
   for(cnti = 0; cnti < Nx; cnti ++) {
      for(cntj = 0; cntj < Ny; cntj ++) {
         psi[cnti][cntj] /= *norm;
      }
   }
   
   return;
}

/**
 *    Calculation of the chemical potential and energy.
 *    mu    - chemical potential
 *    en    - energy
 *    psi   - array with the wave function values
 *    dpsix - temporary array
 *    dpsiy - temporary array
 *    tmpxi - temporary array
 *    tmpyi - temporary array
 *    tmpxj - temporary array
 *    tmpyj - temporary array
 */
void calcmuen(double *mu, double *en, double complex **psi, double **dpsix, double **dpsiy, double *tmpxi, double *tmpyi, double *tmpxj, double *tmpyj) {
   long cnti, cntj;
   double psi2, psi2lin, dpsi2;
   
   for(cntj = 0; cntj < Ny; cntj ++) {
      for(cnti = 0; cnti < Nx; cnti ++) {
         tmpxi[cnti] = cabs(psi[cnti][cntj]);
      }
      diff(dx, tmpxi, tmpxj, Nx);
      for(cnti = 0; cnti < Nx; cnti ++) {
         dpsix[cnti][cntj] = tmpxj[cnti];
      }
   }
   
   for(cnti = 0; cnti < Nx; cnti ++) {
      for(cntj = 0; cntj < Ny; cntj ++) {
         tmpyi[cntj] = cabs(psi[cnti][cntj]);
      }
      diff(dy, tmpyi, tmpyj, Ny);
      for(cntj = 0; cntj < Ny; cntj ++) {
         dpsiy[cnti][cntj] = tmpyj[cntj];
      }
   }
   
   for(cnti = 0; cnti < Nx; cnti ++) {
      for(cntj = 0; cntj < Ny; cntj ++) {
         psi2 = cabs(psi[cnti][cntj]);
         psi2 *= psi2;
         psi2lin = psi2 * G;
         dpsi2 = dpsix[cnti][cntj] * dpsix[cnti][cntj] + 
                 dpsiy[cnti][cntj] * dpsiy[cnti][cntj];
         tmpyi[cntj] = (pot[cnti][cntj] + psi2lin) * psi2 + dpsi2;
         tmpyj[cntj] = (pot[cnti][cntj] + 0.5 * psi2lin) * psi2 + dpsi2;
      }
      tmpxi[cnti] = simpint(dy, tmpyi, Ny);
      tmpxj[cnti] = simpint(dy, tmpyj, Ny);
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
 *    tmpy - temporary array
 */
void calcrms(double *rms, double complex **psi, double *tmpx, double *tmpy) {
   long cnti, cntj;
   double psi2;
   
   for(cnti = 0; cnti < Nx; cnti ++) {
      for(cntj = 0; cntj < Ny; cntj ++) {
         psi2 = cabs(psi[cnti][cntj]);
         psi2 *= psi2;
         tmpy[cntj] = (x2[cnti] + y2[cntj]) * psi2;
      }
      tmpx[cnti] = simpint(dy, tmpy, Ny);
   }
   *rms = sqrt(simpint(dx, tmpx, Nx));
   
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
   
   for(cnti = 0; cnti < Nx; cnti ++) {
      for(cntj = 0; cntj < Ny; cntj ++) {
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
 *    Time propagation with respect to H2 (x-part of the Laplacian).
 *    psi   - array with the wave function values
 *    cbeta - Crank-Nicolson scheme coefficients
 */
void calclux(double complex **psi, double complex *cbeta) {
   long cnti, cntj;
   double complex c;
   
   for(cntj = 0; cntj < Ny; cntj ++) {
      cbeta[Nx - 2] = psi[Nx - 1][cntj];
      for (cnti = Nx - 2; cnti > 0; cnti --) {
         c = - Ax * psi[cnti + 1][cntj] + Ax0r * psi[cnti][cntj] - Ax * psi[cnti - 1][cntj];
         cbeta[cnti - 1] =  cgammax[cnti] * (Ax * cbeta[cnti] - c);
      }
      psi[0][cntj] = 0.;
      for (cnti = 0; cnti < Nx - 2; cnti ++) {
         psi[cnti + 1][cntj] = calphax[cnti] * psi[cnti][cntj] + cbeta[cnti];
      }
      psi[Nx - 1][cntj] = 0.;
   }
   
   return;
}

/**
 *    Time propagation with respect to H3 (y-part of the Laplacian).
 *    psi   - array with the wave function values
 *    cbeta - Crank-Nicolson scheme coefficients
 */
void calcluy(double complex **psi, double complex *cbeta) {
   long cnti, cntj;
   double complex c;
   
   for(cnti = 0; cnti < Nx; cnti ++) {
      cbeta[Ny - 2] = psi[cnti][Ny - 1];
      for (cntj = Ny - 2; cntj > 0; cntj --) {
         c = - Ay * psi[cnti][cntj + 1] + Ay0r * psi[cnti][cntj] - Ay * psi[cnti][cntj - 1];
         cbeta[cntj - 1] =  cgammay[cntj] * (Ay * cbeta[cntj] - c);
      }
      psi[cnti][0] = 0.;
      for (cntj = 0; cntj < Ny - 2; cntj ++) {
         psi[cnti][cntj + 1] = calphay[cntj] * psi[cnti][cntj] + cbeta[cntj];
      }
      psi[cnti][Ny - 1] = 0.;
   }
   
   return;
}
