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
 * 
 *    This program solves the time-dependent Gross–Pitaevskii nonlinear partial
 *    differential equation in one space dimensional trap using the real-time 
 *    propagation. The Gross–Pitaevskii equation describes the properties of a 
 *    dilute trapped Bose–Einstein condensate. The equation is solved using the 
 *    split-step Crank–Nicolson method by discretizing space and time. The
 *    discretized equation is then propagated in real time over small time steps.
 *    
 *    Description of variables used in the code:
 *    
 *    opt     - decides which rescaling of GP equation will be used
 *    par     - parameter for rescaling of GP equation
 *    psi     - array with the wave function values
 *    pot     - array with the values of the potential
 *    g_1d      - final nonlinearity
 *    Gpar    - coefficient that multiplies nonlinear term in non-stationary
 *              problem during final Nrun iterations
 *    norm    - wave function norm
 *    rms     - root mean square radius
 *    mu      - chemical potential
 *    en      - energy
 *    Nx      - number of discretization points in the x-direction
 *    x       - array with the space mesh values in the x-direction
 *    dx      - spatial discretization step in the x-direction
 *    dt      - time discretization step
 *    vgamma  - Gamma coefficient of anisotropy of the trap (omega_x / omega)
 *    Nstp    - number of subsequent iterations with fixed nonlinearity g_1d
 *    Npas    - number of subsequent iterations with the fixed nonlinearity g_1d
 *    Nrun    - number of final iterations with the fixed nonlinearity g_1d
 *    output  - output file with the summary of final values of all physical
 *              quantities
 *    initout - output file with the initial wave function
 *    Nstpout - output file with the wave function obtained after the first
 *              Nstp iterations
 *    Npasout - output file with the wave function obtained after the
 *              subsequent Npas iterations, with the fixed nonlinearity g_1d
 *    Nrunout - output file with the final wave function obtained after the
 *              final Nrun iterations
 *    dynaout  - output file with the time dependence of RMS during Npas
 *              and the final Nrun iterations
 *    outstpx - discretization step in the x-direction used to save wave
 *              functions
 *    outstpt - time discretization step used to save RMS of the wave function
 */

#include "real1d.h"

int main(int argc, char **argv) {
   FILE *out;
   FILE *file;
   FILE *dyna;
   char filename[MAX_FILENAME_SIZE];
   long cnti;
   double norm, mu, en;
   double *rms;
   double complex *psi;
   double complex *cbeta;
   double *dpsix;
   double *tmpxi, *tmpxj;
   double psi2;
   double *abc;
      
   pi = 3.14159265358979;
   
   time_t clock_beg, clock_end; 
   clock_beg = time(NULL);

   if((argc != 3) || (strcmp(*(argv + 1), "-i") != 0)) {
      fprintf(stderr, "Usage: %s -i <input parameter file> \n", *argv);
      exit(EXIT_FAILURE);
   }

   if(! cfg_init(argv[2])) {
      fprintf(stderr, "Wrong input parameter file.\n");
      exit(EXIT_FAILURE);
   }

   readpar();

   rms = alloc_double_vector(RMS_ARRAY_SIZE);

   x = alloc_double_vector(Nx);
   x2 = alloc_double_vector(Nx);
   abc = alloc_double_vector(Nx);

   pot = alloc_double_vector(Nx);
   psi = alloc_complex_vector(Nx);

   dpsix = alloc_double_vector(Nx);

   calphax = alloc_complex_vector(Nx - 1);
   cbeta =  alloc_complex_vector(Nx - 1);
   cgammax = alloc_complex_vector(Nx - 1);

   tmpxi = alloc_double_vector(Nx);
   tmpxj = alloc_double_vector(Nx);

   if(output != NULL) {
      sprintf(filename, "%s.txt", output);
      out = fopen(filename, "w");
   }
   else out = stdout;

   if (opt == 1) par = 1.;
   else if (opt == 2) par = 2.;
   else{
      fprintf(stderr, "OPTION is not well defined in the configuration file\n");
      exit(EXIT_FAILURE);
   }

   Nx2 = Nx / 2;
   dx2 = dx * dx;

   fprintf(out, " Real time propagation 1d,   OPTION = %d\n\n", opt);
   fprintf(out, "  Number of Atoms N = %li, Unit of length AHO = %.8f m\n", Na,aho);
   fprintf(out, "  Scattering length a = %.2f*a0\n", as);
   fprintf(out, "  Nonlinearity G_3D = %.7f\n", g_3d);
   fprintf(out, "  Nonlinearity G_1D = %.7f\n", g_1d);
   fprintf(out, "  Parameter of trap: GAMMA = %.2f\n", vgamma);
   fprintf(out, "  Radial Gaussian Width = %.2f\n\n", drho);
   fprintf(out, " # Space Stp:  N = %li\n", Nx);
   fprintf(out, "  Space Step: DX = %.6f\n", dx);
   fprintf(out, " # Time Stp : NSTP = %li, NPAS = %li, NRUN = %li\n", Nstp, Npas, Nrun);
   fprintf(out, "   Time Step:   DT = %.6f\n\n",  dt);
   fprintf(out, " * Change for dynamics: GPAR = %.3f *\n\n", gpar);   
   fprintf(out, "                  --------------------------------------------------------\n");
   fprintf(out, "                    Norm      Chem        Ener/N      <x>      |Psi(0)|^2\n");
   fprintf(out, "                  --------------------------------------------------------\n");
   fflush(out);

   if(Nstp == 0) { 
     g = par * g_1d;
   }
   else {
     g = 0.;
   }

   initpsi(psi, abc);
   initpot();
   gencoef();   
   calcnorm(&norm, psi, tmpxi);
   calcmuen(&mu, &en, psi, dpsix, tmpxi, tmpxj);
   calcrms(rms, psi, tmpxi);
   psi2 = cabs(psi[Nx2]);
   psi2 *= psi2;
   fprintf(out, "Initial : %15.4f %11.4f %11.4f %10.3f %10.3f\n", norm, mu / par, en / par, *rms, psi2);
   fflush(out);
   if(initout != NULL) {
      sprintf(filename, "%s.txt", initout);
      file = fopen(filename, "w");
      outpsi2x(psi, file);
      fclose(file);
   }
   
   if(Nstp != 0) {
     double g_stp = par * g_1d / (double) Nstp;
     for(cnti = 0; cnti < Nstp; cnti ++) {
        g += g_stp;
        calcnu(psi);
        calclux(psi, cbeta);
     }
    calcnorm(&norm, psi, tmpxi);
    calcmuen(&mu, &en, psi, dpsix, tmpxi, tmpxj);
    calcrms(rms, psi, tmpxi);
    psi2 = cabs(psi[Nx2]);
    psi2 *= psi2;
    fprintf(out, "After NSTP iter.:%8.4f %11.4f %11.4f %10.3f %10.3f\n", norm, mu / par, en / par, *rms, psi2);
    fflush(out);
    if(Nstpout != NULL) {
      sprintf(filename, "%s.txt", Nstpout);
      file = fopen(filename, "w");
      outpsi2x(psi, file);
      fclose(file);
    }
   }
   
   if(dynaout != NULL) {
      sprintf(filename, "%s.txt", dynaout);
      dyna = fopen(filename, "w");
   }
   else dyna = NULL;

   for(cnti = 0; cnti < Npas; cnti ++) {
      calcnu(psi);
      calclux(psi, cbeta);
      if((dynaout != NULL) && (cnti % outstpt == 0)) {
         calcrms(rms, psi, tmpxi);
         fprintf(dyna, "%8le   %8le\n", cnti * dt * par, *rms);
         fflush(dyna);
      }
   }
   if(Npas != 0){    
     calcnorm(&norm, psi, tmpxi);
     calcmuen(&mu, &en, psi, dpsix, tmpxi, tmpxj);
     calcrms(rms, psi, tmpxi);
     psi2 = cabs(psi[Nx2]);
     psi2 *= psi2;
     fprintf(out, "After NPAS iter.:%8.4f %11.4f %11.4f %10.3f %10.3f\n", norm, mu / par, en / par, *rms, psi2);
     fflush(out);
     if(Npasout != NULL) {
        sprintf(filename, "%s.txt", Npasout);
        file = fopen(filename, "w");
        outpsi2x(psi, file);
        fclose(file);
     }
   }

   g *= gpar;

   for(cnti = 0; cnti < Nrun; cnti ++) {
      calcnu(psi);
      calclux(psi, cbeta);
      if((dynaout != NULL) && (cnti % outstpt == 0)) {
         calcrms(rms, psi, tmpxi);
         fprintf(dyna, "%8le   %8le\n", (cnti * dt  + dt * Npas) * par, *rms);
         fflush(dyna);
      }
   }
   if(dynaout != NULL) fclose(dyna);
   if(Nrun != 0){   
     calcnorm(&norm, psi, tmpxi);
     calcmuen(&mu, &en, psi, dpsix, tmpxi, tmpxj);
     calcrms(rms, psi, tmpxi);
     psi2 = cabs(psi[Nx2]);
     psi2 *= psi2;
     fprintf(out, "After NRUN iter.:%8.4f %11.4f %11.4f %10.3f %10.3f\n", norm, mu / par, en / par, *rms, psi2);
     fflush(out);
     if(Nrunout != NULL) {
        sprintf(filename, "%s.txt", Nrunout);
        file = fopen(filename, "w");
        outpsi2x(psi, file);
        fclose(file);
     }
   }
   fprintf(out, "                  --------------------------------------------------------\n");

   free_double_vector(rms);

   free_double_vector(x);

   free_double_vector(x2);

   free_double_vector(pot);
   free_double_vector(abc);
   free_complex_vector(psi);

   free_double_vector(dpsix);

   free_complex_vector(calphax);
   free_complex_vector(cbeta);
   free_complex_vector(cgammax);

   free_double_vector(tmpxi);
   free_double_vector(tmpxj);
   
   clock_end = time(NULL);
   //fprintf(out, "Clock start : %s", ctime(&clock_beg));
   //fprintf(out, "Clock end   : %s", ctime(&clock_end));
   double wall_time = difftime(clock_end, clock_beg);   
   double cpu_time = clock() / (double) CLOCKS_PER_SEC;
   fprintf(out, " Clock Time: %.f seconds\n", wall_time);
   fprintf(out, " CPU Time: %.f seconds\n", cpu_time);  
   
   if(output != NULL) fclose(out);

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

   if((cfg_tmp = cfg_read("DRHO")) == NULL) {
      fprintf(stderr, "DRHO is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   drho = atof(cfg_tmp); 
   drho2 = drho * drho;

   if((cfg_tmp = cfg_read("G_1D")) == NULL) {
      if((cfg_tmp = cfg_read("NATOMS")) == NULL) {
	 fprintf(stderr, "NATOMS is not defined in the configuration file.\n");
	 exit(EXIT_FAILURE);
      }
      Na = atol(cfg_tmp);
   
      if((cfg_tmp = cfg_read("AHO")) == NULL) {
         fprintf(stderr, "AHO is not defined in the configuration file.\n");
         exit(EXIT_FAILURE);
      }
      aho = atof(cfg_tmp);

      if((cfg_tmp = cfg_read("AS")) == NULL) {
         fprintf(stderr, "AS is not defined in the configuration file.\n");
         exit(EXIT_FAILURE);
      }
      as = atof(cfg_tmp);

      g_3d = 4. * pi * as * Na * BOHR_RADIUS / aho;
      g_1d = g_3d / (2. * pi * drho2);      
   } else {
     g_1d = atof(cfg_tmp); 
     g_3d = g_1d * (2. * pi * drho2);
   }

   if((cfg_tmp = cfg_read("GPAR")) == NULL) {
      fprintf(stderr, "GPAR is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   gpar = atof(cfg_tmp);

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

   if((cfg_tmp = cfg_read("GAMMA")) == NULL) {
      fprintf(stderr, "GAMMA is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   vgamma = atof(cfg_tmp);

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
   dynaout = cfg_read("DYNAOUT");
   initout = cfg_read("INITOUT");
   Nstpout = cfg_read("NSTPOUT");
   Npasout = cfg_read("NPASOUT");
   Nrunout = cfg_read("NRUNOUT");

   if((initout != NULL) || (Nstpout != NULL) || (Npasout != NULL) || (Nrunout != NULL)) {
      if((cfg_tmp = cfg_read("OUTSTPX")) == NULL) {
         fprintf(stderr, "OUTSTPX is not defined in the configuration file.\n");
         exit(EXIT_FAILURE);
      }
      outstpx = atol(cfg_tmp);
   }

   if(dynaout != NULL) {
      if((cfg_tmp = cfg_read("OUTSTPT")) == NULL) {
         fprintf(stderr, "OUTSTPT is not defined in the configuration file.\n");
         exit(EXIT_FAILURE);
      }
      outstpt = atol(cfg_tmp);
   }

   return;
}

/**
 *    Initialization of the space mesh and the initial wave function.
 *    psi - array with the wave function values
 */
void initpsi(double complex *psi, double *abc) {
   long cnti;
   double cpsi;
   double tmp;
   FILE *file;

   cpsi = sqrt( sqrt(pi  / vgamma));

   for(cnti = 0; cnti < Nx; cnti ++) {
      x[cnti] = (cnti - Nx2) * dx;
      x2[cnti] = x[cnti] * x[cnti];
   }
   
   if(Nstp!=0  ) {
      for(cnti = 0; cnti < Nx; cnti ++) {
	  tmp = exp(- 0.5 * vgamma * x2[cnti]);
	  psi[cnti] = tmp / cpsi;
      }
   }
   
   if(Nstp==0  ) {
     if((file = fopen("imag1d-den.txt", "r"))==NULL) {       /* open a text file for reading */
       printf("Run the program using the input file to read. i.g.: imag1d-den.txt\n");
       exit(1);					/*couldn't open the requested file!*/
     }
     file = fopen("imag1d-den.txt", "r");
     for(cnti = 0; cnti < Nx; cnti ++) {
	if(fscanf(file,"%lf %lf\n",&x[cnti],&abc[cnti])) {
//            printf("%8le %8le\n",x[cnti], abc[cnti]);
	}
     }
     fclose(file);
     for(cnti = 0; cnti < Nx; cnti ++) {
       psi[cnti]=sqrt(abc[cnti]);
     }
   }

   return;
}

/**
 *    Initialization of the potential.
 */
void initpot() {
   long cnti;
   double vgamma2;

   vgamma2 = vgamma * vgamma;

   for(cnti = 0; cnti < Nx; cnti ++) {
      pot[cnti] = (vgamma2 * x2[cnti]);
   }

   return;
}

/**
 *    Crank-Nicolson scheme coefficients generation.
 */
void gencoef(void) {
   long cnti;

   Ax0 = 1. + I * dt / dx2;

   Ax0r = 1. - I * dt / dx2;

   Ax = - 0.5 * I * dt / dx2;

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
 *    tmpy - temporary array
 */
void calcnorm(double *norm, double complex *psi, double *tmpx) {
   long cnti;
   //double tmp;

   for(cnti = 0; cnti < Nx; cnti ++) {
      tmpx[cnti] = cabs(psi[cnti]);
      tmpx[cnti] *= tmpx[cnti];
   }

   *norm = sqrt(simpint(dx, tmpx, Nx));
   //tmp = 1. / *norm;

   //for(cnti = 0; cnti < Nx; cnti ++) {
   //   psi[cnti] *= tmp;
   //}

   return;
}

/**
 *    Calculation of the chemical potential and energy.
 *    mu        - chemical potential
 *    en        - energy
 *    psi       - array with the wave function values
 *    dpsix     - temporary array
 *    dpsiy     - temporary array
 *    tmpxi     - temporary array
 *    tmpyi     - temporary array
 *    tmpzi     - temporary array
 *    tmpxj     - temporary array
 *    tmpyj     - temporary array
 */
void calcmuen(double *mu, double *en, double complex *psi, double *dpsix, double *tmpxi, double *tmpxj) {
   long cnti;
   double psi2, psi2lin, dpsi2;

   for(cnti = 0; cnti < Nx; cnti ++) {
      tmpxi[cnti] = cabs(psi[cnti]);
   }  

   diff(dx, tmpxi, dpsix, Nx);

   for(cnti = 0; cnti < Nx; cnti ++) {
      psi2 = cabs(psi[cnti]);
      psi2 *= psi2;
      psi2lin = psi2 * g;
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
 *    tmpy - temporary array
 */
void calcrms(double *rms, double complex *psi, double *tmpx) {
   long cnti;
   double psi2;

   for(cnti = 0; cnti < Nx; cnti ++) {
      psi2 = cabs(psi[cnti]);
      psi2 *= psi2;
      tmpx[cnti] = x2[cnti] * psi2;
   }
   *rms = sqrt(simpint(dx, tmpx, Nx));

   return;
}

/**
 *    Time propagation with respect to H1 (part of the Hamiltonian without 
 *    spatial derivatives).
 *    psi       - array with the wave function values
 */
void calcnu(double complex *psi) {
   long cnti;
   double psi2, psi2lin, tmp;

   #pragma omp parallel for private(cnti, psi2, psi2lin, tmp)    
   for(cnti = 0; cnti < Nx; cnti ++) {
      psi2 = cabs(psi[cnti]);
      psi2 *= psi2;
      psi2lin = psi2 * g;
      tmp = dt * (pot[cnti] + psi2lin);
      psi[cnti] *= cexp(- I * tmp);
   }

   return;
}

/**
 *    Time propagation with respect to H2 (x-part of the Laplacian).
 *    psi   - array with the wave function values
 *    cbeta - Crank-Nicolson scheme coefficients
 */
void calclux(double complex *psi, double complex *cbeta) {
   long cnti;
   double complex c;

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

void outpsi2x(double complex *psi, FILE *file) {
   long cnti;

   for(cnti = 0; cnti < Nx; cnti += outstpx) {
      fprintf(file, "%8le %8le\n", x[cnti], cabs(psi[cnti]) * cabs(psi[cnti]));
   }
   fflush(file);

   return;
}
