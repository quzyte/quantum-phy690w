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
 *    This program solves the time-independent Gross–Pitaevskii nonlinear
 *    partial differential equation in two space dimensions in a
 *    circularly-symmetric trap using the imaginary-time propagation. The
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
 *    g_cir   - final nonlinearity
 *    norm    - wave function norm
 *    rms     - root mean square radius
 *    mu      - chemical potential
 *    en      - energy
 *    Nr      - number of discretization points in the r-direction
 *    r       - array with the space mesh values in the r-direction
 *    dr      - spatial discretization step in the r-direction
 *    dt      - time discretization step
 *    Nstp    - number of initial iterations to introduce the nonlinearity g_cir
 *    Npas    - number of subsequent iterations with the fixed nonlinearity g_cir
 *    Nrun    - number of final iterations with the fixed nonlinearity g_cir
 *    output  - output file with the summary of final values of all physical 
 *              quantities
 *    initout - output file with the initial wave function
 *    Npasout - output file with the wave function obtained after the 
 *              subsequent Npas iterations, with the fixed nonlinearity G_CIR
 *    Nrunout - output file with the final wave function obtained after the 
 *              final Nrun iterations
 *    outstpr - discretization step in the r-direction used to save wave 
 *              functions
 */

#include "imagcir.h"

int main(int argc, char **argv) {
   FILE *out;
   FILE *file;
   char filename[MAX_FILENAME_SIZE];    
   long cnti;
   double norm, rms, mu, en;
   double *psi;
   double *cbeta;
   double *dpsir;
   double *tmpri, *tmprj;
   
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
   
   r = alloc_double_vector(Nr);
   r2 = alloc_double_vector(Nr);
   r3 = alloc_double_vector(Nr);

   pot = alloc_double_vector(Nr);
   psi = alloc_double_vector(Nr);
   
   dpsir = alloc_double_vector(Nr);
   
   calphar = alloc_double_vector(Nr - 1);
   cbeta =  alloc_double_vector(Nr - 1);
   cgammar = alloc_double_vector(Nr - 1);
   
   tmpri = alloc_double_vector(Nr);
   tmprj = alloc_double_vector(Nr);

   if(output != NULL) {
      sprintf(filename, "%s.txt", output);
      out = fopen(filename, "w");
   }
   else out = stdout;
   
   fprintf(out, " Imaginary time propagation circularly-symmetric trap,   OPTION = %d\n\n", opt);
   fprintf(out, "  Number of Atoms N = %li, Unit of length AHO = %.8f m\n", Na, aho);
   fprintf(out, "  Scattering length a = %.2f*a0\n", as);
   fprintf(out, "  Nonlinearity G_3D = %.7f\n", g_3d);
   fprintf(out, "  Nonlinearity G_2D = %.7f\n", g_cir);
   fprintf(out, "  Axial trap parameter = %.2f\n\n", d_z);
   fprintf(out, " # Space Stp: NR = %li\n", Nr);
   fprintf(out, "  Space Step: DR = %.4f\n", dr);
   fprintf(out, " # Time Stp : NSTP = %li, NPAS = %li, NRUN = %li\n", Nstp, Npas, Nrun);
   fprintf(out, "   Time Step:   DT = %.6f\n\n",  dt);
   fprintf(out, "                  --------------------------------------------------------\n");
   fprintf(out, "                    Norm       Chem       Ener/N      <rms>    |Psi(0)|^2\n");
   fprintf(out, "                  --------------------------------------------------------\n");
   fflush(out);   
   
   G = 0.;
   
   init(psi);
   gencoef();
   calcnorm(&norm, psi, tmpri);
   calcmuen(&mu, &en, psi, dpsir, tmpri, tmprj);
   calcrms(&rms, psi, tmpri);
   fprintf(out, "Initial : %15.4f %11.6f %11.6f %10.5f %10.5f\n", norm, mu / par, en / par, rms, psi[0] * psi[0]);
   fflush(out);
   
   if(initout != NULL) {
      sprintf(filename, "%s.txt", initout);
      file = fopen(filename, "w");
      for(cnti = 0; cnti < Nr; cnti += outstpr)
         fprintf(file, "%8le %8le\n", r[cnti], psi[cnti] * psi[cnti]);
      fclose(file);
   }
   
   if(Nstp != 0) {
     double g_stp = par * g_cir / (double) Nstp;
     G = 0.;
     for(cnti = 0; cnti < Nstp; cnti ++) {
        G += g_stp;
	calcnu(psi);
	calclur(psi, cbeta);
	calcnorm(&norm, psi, tmpri);
     }
     calcmuen(&mu, &en, psi, dpsir, tmpri, tmprj);
     calcrms(&rms, psi, tmpri);
     fprintf(out, "After NSTP iter.:%8.4f %11.6f %11.6f %10.5f %10.5f\n", norm, mu / par, en / par, rms, psi[0] * psi[0]);
     fflush(out);
     
     if(Nstpout != NULL) {
       sprintf(filename, "%s.txt", Nstpout);
       file = fopen(filename, "w");
      for(cnti = 0; cnti < Nr; cnti += outstpr)
         fprintf(file, "%8le %8le\n", r[cnti], psi[cnti] * psi[cnti]);
       fclose(file);
     }
   }
   else {
     G = par * g_cir;
   }
   
   for(cnti = 0; cnti < Npas; cnti ++) {
      calcnu(psi);
      calclur(psi, cbeta);
      calcnorm(&norm, psi, tmpri);
   }
   if(Npas != 0){    
    calcmuen(&mu, &en, psi, dpsir, tmpri, tmprj);
    calcrms(&rms, psi, tmpri);
    fprintf(out, "After NPAS iter.:%8.4f %11.6f %11.6f %10.5f %10.5f\n", norm, mu / par, en / par, rms, psi[0] * psi[0]);
   fflush(out);
     
    if(Npasout != NULL) {
      sprintf(filename, "%s.txt", Npasout);
      file = fopen(filename, "w");
      for(cnti = 0; cnti < Nr; cnti += outstpr)
         fprintf(file, "%8le %8le\n", r[cnti], psi[cnti] * psi[cnti]);
      fclose(file);
    }
   }
   
   for(cnti = 0; cnti < Nrun; cnti ++) {
      calcnu(psi);
      calclur(psi, cbeta);
      calcnorm(&norm, psi, tmpri);
   }
   if(Nrun != 0){    
    calcmuen(&mu, &en, psi, dpsir, tmpri, tmprj);
    calcrms(&rms, psi, tmpri);
    fprintf(out, "After NRUN iter.:%8.4f %11.6f %11.6f %10.5f %10.5f\n", norm, mu / par, en / par, rms, psi[0] * psi[0]);
   fflush(out);
   
    if(Nrunout != NULL) {
      sprintf(filename, "%s.txt", Nrunout);
      file = fopen(filename, "w");
       for(cnti = 0; cnti < Nr; cnti += outstpr)
          fprintf(file, "%8le %8le\n", r[cnti], psi[cnti] * psi[cnti]);
       fclose(file);
    }
   }
   
   fprintf(out, "                  --------------------------------------------------------\n\n");
   
   free_double_vector(r);
   free_double_vector(r2);
   free_double_vector(r3);

   free_double_vector(pot);
   free_double_vector(psi);

   free_double_vector(dpsir);

   free_double_vector(calphar);
   free_double_vector(cbeta);
   free_double_vector(cgammar);

   free_double_vector(tmpri);
   free_double_vector(tmprj);
   
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
   
   if((cfg_tmp = cfg_read("D_Z")) == NULL) {
      fprintf(stderr, "D_Z is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   d_z = atof(cfg_tmp);    
   
   if((cfg_tmp = cfg_read("G_CIR")) == NULL) {
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

      g_3d = 4. * (4. * atan(1.)) * as * Na * BOHR_RADIUS / aho;     
      g_cir = g_3d / (sqrt(2. * (4. * atan(1.))) * d_z);
   } else {
      g_cir = atof(cfg_tmp);
      g_3d = g_cir * (sqrt(2. * (4. * atan(1.))) * d_z);      
   }

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
   Nstpout = cfg_read("NSTPOUT");   
   Npasout = cfg_read("NPASOUT");
   Nrunout = cfg_read("NRUNOUT");
   
   if((initout != NULL) || (Nstpout != NULL) || (Npasout != NULL) || (Nrunout != NULL)) {
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
   double cpsi1;
   double tmp;
   
   if (opt == 1) par = 1.;
   else if (opt == 2) par = 2.;
   else{
      fprintf(stderr, "OPTION is not well defined in the configuration file\n");
      exit(EXIT_FAILURE);
   }
   
   Nr2 = Nr / 2;
   dr2 = dr * dr;
   
   pi = 4. * atan(1.);
   cpsi1 = sqrt(pi);
   
   for(cnti = 0; cnti < Nr; cnti ++) {
      r[cnti] = cnti * dr;
      r2[cnti] = r[cnti] * r[cnti];
      r3[cnti] = r2[cnti] * r[cnti];

      pot[cnti] = r2[cnti];
      tmp = exp(- 0.5 * r2[cnti]);
      psi[cnti] = tmp / cpsi1;      
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
   
   dAr = 0.25 * dt / dr;

   calphar[0] = 1.;
   cgammar[0] = - 1. / (Ar0 + Ar + dAr / r[1]);
   for(cnti = 1; cnti < Nr - 1; cnti ++) {
      calphar[cnti] = (Ar - dAr / r[cnti]) * cgammar[cnti - 1];
      cgammar[cnti] = - 1. / (Ar0 + (Ar + dAr / r[cnti + 1]) * calphar[cnti]);
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
	  tmpr[cnti] = psi[cnti] * psi[cnti] * r[cnti];
      }

      *norm = sqrt(2. * pi * simpint(dr, tmpr, Nr));
   
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

      for(cnti = 0; cnti < Nr; cnti ++) {
	  psi2 = psi[cnti] * psi[cnti];
	  psi2lin = psi2 * G;
	  dpsi2 = dpsir[cnti] * dpsir[cnti];
	  tmpri[cnti] = r[cnti] * ((pot[cnti] + psi2lin) * psi2 + dpsi2);
	  tmprj[cnti] = r[cnti] * ((pot[cnti] + 0.5 * psi2lin) * psi2 + dpsi2);
      }

   *mu = 2. * pi * simpint(dr, tmpri, Nr);
   *en = 2. * pi * simpint(dr, tmprj, Nr);
   
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
      tmpr[cnti] = r3[cnti] * psi2;
   } 
   *rms = sqrt(2. * pi * simpint(dr, tmpr, Nr));
   
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

   #pragma omp parallel for private(cnti, psi2, psi2lin, tmp)   
   for(cnti = 0; cnti < Nr; cnti ++) {
      psi2 = psi[cnti] * psi[cnti];
      psi2lin = psi2 * G;
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
   
   cbeta[0] = 0.;
   for(cnti = 1; cnti < Nr - 1; cnti ++) {
      c = psi[cnti] - Ar * (psi[cnti + 1] - 2 * psi[cnti] + psi[cnti - 1]) + dAr / r[cnti] * (psi[cnti + 1] - psi[cnti - 1]);
      cbeta[cnti] = cgammar[cnti - 1] * ((Ar + dAr / r[cnti]) * cbeta[cnti - 1] - c);
   }
   psi[Nr - 1] = 0.;
   for (cnti = Nr - 1; cnti > 0; cnti --) {
      psi[cnti - 1] = calphar[cnti - 1] * psi[cnti] + cbeta[cnti - 1];
   }
   
   return;
}
