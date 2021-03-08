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
 *    partial differential equation in three space dimensions in an
 *    axially-symmetric trap using the imaginary-time propagation. The
 *    Gross–Pitaevskii equation describes the properties of a dilute trapped
 *    Bose–Einstein condensate. The equation is solved using the split-step
 *    Crank–Nicolson method by discretizing space and time. The discretized
 *    equation is then propagated in imaginary time over small time steps.
 *    When convergence is achieved, the method has yielded the stationary
 *    solution of the problem.
 *    
 *    Description of variables used in the code:
 *    
 *   opt       - decides which rescaling of GP equation will be used
 *   par       - parameter for rescaling of GP equation
 *   psi       - array with the wave function values
 *   pot       - array with the values of the potential
 *   G0        - final nonlinearity
 *   norm      - wave function norm
 *   rms       - root mean square radius
 *   mu        - chemical potential
 *   en        - energy
 *   Nrho      - number of discretization points in the rho-direction
 *   Nz        - number of discretization points in the z-direction
 *   rho       - array with the space mesh values in the rho-direction
 *   z         - array with the space mesh values in the z-direction
 *   drho      - spatial discretization step in the rho-direction
 *   dz        - spatial discretization step in the z-direction
 *   dt        - time discretization step
 *   vnu       - nu coefficient of anisotropy of the trap
 *   vlambda   - lambda coefficient of anisotropy of the trap
 *   Nstp      - number of initial iterations to introduce the nonlinearity G0
 *   Npas      - number of subsequent iterations with the fixed nonlinearity G0
 *   Nrun      - number of final iterations with the fixed nonlinearity G0
 *   output    - output file with the summary of final values of all physical
 *               quantities
 *   initout   - output file with the initial wave function
 *   Npasout   - output file with the wave function obtained after the
 *               subsequent Npas iterations, with the fixed nonlinearity G0
 *   Nrunout   - output file with the final wave function obtained after the
 *               final Nrun iterations
 *   outstprho - discretization step in the rho-direction used to save wave
 *               functions
 *   outstpz   - discretization step in the z-direction used to save wave
 *               functions
 */

#include "imagaxi.h"

int main(int argc, char **argv) {
   FILE *out;
   FILE *file;
   FILE *filerms;      
   char filename[MAX_FILENAME_SIZE];    
   int nthreads;
   long cnti, cntj;
   double norm, rmsrho, rmsz, rms, mu, en;
   double **psi;
   double **cbeta;
   double **dpsirho, **dpsiz;
   double **tmprhoi, **tmpzi, **tmprhoj, **tmpzj;
   double *tmprho, *tmpz;   
   double psi2;   
   
   time_t clock_beg, clock_end; 
   clock_beg = time(NULL);      
   
   pi = 4. * atan(1.);
   
   if((argc != 3) || (strcmp(*(argv + 1), "-i") != 0)) {
      fprintf(stderr, "Usage: %s -i <input parameter file> \n", *argv);
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
   psi = alloc_double_matrix(Nrho, Nz);
   
   dpsirho = alloc_double_matrix(Nrho, Nz);
   dpsiz = alloc_double_matrix(Nrho, Nz);
   
   calpharho = alloc_double_vector(Nrho - 1);
   calphaz = alloc_double_vector(Nz - 1);
   cbeta =  alloc_double_matrix(nthreads, MAX(Nrho, Nz) - 1);
   cgammarho = alloc_double_vector(Nrho - 1);
   cgammaz = alloc_double_vector(Nz - 1);
   
   tmprhoi = alloc_double_matrix(nthreads, Nrho);
   tmpzi = alloc_double_matrix(nthreads, Nz);
   tmprhoj = alloc_double_matrix(nthreads, Nrho);
   tmpzj = alloc_double_matrix(nthreads, Nz);
   
   tmprho = alloc_double_vector(Nrho);
   tmpz = alloc_double_vector(Nz);
   
   if(output != NULL) {
      sprintf(filename, "%s.txt", output);
      out = fopen(filename, "w");
   }
   else out = stdout;
   
   if(rmsout != NULL) {
      sprintf(filename, "%s.txt", rmsout);
      filerms = fopen(filename, "w");
   }   
   else filerms = NULL;      
   
   fprintf(out, " Imaginary time propagation axially-symmetric trap,   OPTION = %d\n\n", opt);
   fprintf(out, "  Number of Atoms N = %li, Unit of length AHO = %.8f m\n", Na, aho);
   fprintf(out, "  Scattering length a = %.2f*a0\n", as);
   fprintf(out, "  Nonlinearity G_3D = %.7f\n", G0);
   fprintf(out, "  Parameters of trap: NU = %.2f, LAMBDA = %.2f\n", vnu, vlambda);
   fprintf(out, " # Space Stp: NRHO = %li, NZ = %li\n", Nrho, Nz);
   fprintf(out, "  Space Step: DRHO = %.4f, DZ = %.4f\n", drho, dz);
   fprintf(out, " # Time Stp : NSTP = %li, NPAS = %li, NRUN = %li\n", Nstp, Npas, Nrun);
   fprintf(out, "   Time Step:   DT = %.6f\n\n",  dt);
   fprintf(out, "                  --------------------------------------------------------\n");
   fprintf(out, "                    Norm      Chem        Ener/N     <r>     |Psi(0,0)|^2\n");
   fprintf(out, "                  --------------------------------------------------------\n");
   fflush(out);   
   
   init(psi);
   gencoef();
   calcnorm(&norm, psi, tmprhoi, tmpzi);
   calcmuen(&mu, &en, psi, dpsirho, dpsiz, tmprhoi, tmpzi, tmprhoj, tmpzj);
   calcrms(&rmsrho, &rmsz, &rms, psi, tmprhoi, tmpzi, tmprhoj, tmpzj);
   psi2 = psi[0][Nz2 + 1] * psi[0][Nz2 + 1];
   fprintf(out, "Initial : %15.4f %11.5f %11.5f %10.5f %10.5f\n", norm, mu / par, en / par, rms, psi2);
   fflush(out);     
   
   if(rmsout != NULL) {
      fprintf(filerms, " Imaginary time propagation axially-symmetric trap,   OPTION = %d\n\n", opt);
      fprintf(filerms, "                  --------------------------------------------------------\n");
      fprintf(filerms, "Values of rms size:        <rho>             <z>             <r>\n");
      fprintf(filerms, "                  --------------------------------------------------------\n");      
      fprintf(filerms, "           Initial:%14.5f   %13.5f   %13.5f\n", rmsrho, rmsz, rms);
      fflush(filerms);
   }    
   
   if(initout != NULL) {
      sprintf(filename, "%s.txt", initout);
      file = fopen(filename, "w");
      outpsi2xy(psi, file);
      fclose(file);

      sprintf(filename, "%s1d_rho.txt", initout);
      file = fopen(filename, "w");
      outdenrho(psi, tmpz, file);
      fclose(file);

      sprintf(filename, "%s1d_z.txt", initout);
      file = fopen(filename, "w");
      outdenz(psi, tmprho, file);
      fclose(file);   
   }  
   
   if(Nstp != 0) {
     double g_stp = par * G0 / (double) Nstp;
     G = 0.;
     for(cnti = 0; cnti < Nstp; cnti ++) {
        G += g_stp;
	calcnu(psi);
	calclurho(psi, cbeta);
	calcluz(psi, cbeta);
	calcnorm(&norm, psi, tmprhoi, tmpzi);
     }

     calcmuen(&mu, &en, psi, dpsirho, dpsiz, tmprhoi, tmpzi, tmprhoj, tmpzj);
     calcrms(&rmsrho, &rmsz, &rms, psi, tmprhoi, tmpzi, tmprhoj, tmpzj);
     psi2 = psi[0][Nz2 + 1] * psi[0][Nz2 + 1];
     fprintf(out, "After NSTP iter.:%8.4f %11.5f %11.5f %10.5f %10.5f\n", norm, mu / par, en / par, rms, psi2);
     fflush(out);    
     if(rmsout != NULL) {
       fprintf(filerms, "  After NSTP iter.:%14.5f   %13.5f   %13.5f\n", rmsrho, rmsz, rms);
       fflush(filerms);
     }     

     if(Nstpout != NULL) {
       sprintf(filename, "%s.txt", Nstpout);
       file = fopen(filename, "w");
       outpsi2xy(psi, file);
       fclose(file);

       sprintf(filename, "%s1d_rho.txt", Nstpout);
       file = fopen(filename, "w");
       outdenrho(psi, tmpz, file);
       fclose(file);     

       sprintf(filename, "%s1d_z.txt", Nstpout);
       file = fopen(filename, "w");
       outdenz(psi, tmprho, file);
       fclose(file);        
     }     
   }
   else {   
     G = par * G0;
   }
   
   for(cnti = 0; cnti < Npas; cnti ++) {
      calcnu(psi);
      calclurho(psi, cbeta);
      calcluz(psi, cbeta);
      calcnorm(&norm, psi, tmprhoi, tmpzi);
   }
   if(Npas != 0){   
      calcmuen(&mu, &en, psi, dpsirho, dpsiz, tmprhoi, tmpzi, tmprhoj, tmpzj);
      calcrms(&rmsrho, &rmsz, &rms, psi, tmprhoi, tmpzi, tmprhoj, tmpzj);
      psi2 = psi[0][Nz2 + 1] * psi[0][Nz2 + 1];
      fprintf(out, "After NPAS iter.:%8.4f %11.5f %11.5f %10.5f %10.5f\n", norm, mu / par, en / par, rms, psi2);    
      fflush(out);      
      if(rmsout != NULL) {
	fprintf(filerms, "  After NPAS iter.:%14.5f   %13.5f   %13.5f\n", rmsrho, rmsz, rms);
	fflush(filerms);
      }       
     
      if(Npasout != NULL) {
	sprintf(filename, "%s.txt", Npasout);
	file = fopen(filename, "w");
	outpsi2xy(psi, file);
	fclose(file);

	sprintf(filename, "%s1d_rho.txt", Npasout);
	file = fopen(filename, "w");
	outdenrho(psi, tmpz, file);
	fclose(file);     

	sprintf(filename, "%s1d_z.txt", Npasout);
	file = fopen(filename, "w");
	outdenz(psi, tmprho, file);
	fclose(file);     
      }
   }
   
   for(cnti = 0; cnti < Nrun; cnti ++) {
      calcnu(psi);
      calclurho(psi, cbeta);
      calcluz(psi, cbeta);
      calcnorm(&norm, psi, tmprhoi, tmpzi);
   }
   if(Nrun != 0){ 
      calcmuen(&mu, &en, psi, dpsirho, dpsiz, tmprhoi, tmpzi, tmprhoj, tmpzj);
      calcrms(&rmsrho, &rmsz, &rms, psi, tmprhoi, tmpzi, tmprhoj, tmpzj);
      psi2 = psi[0][Nz2 + 1] * psi[0][Nz2 + 1];
      fprintf(out, "After NRUN iter.:%8.4f %11.5f %11.5f %10.5f %10.5f\n", norm, mu / par, en / par, rms, psi2);         
      fflush(out);
      
      if(rmsout != NULL) {
	fprintf(filerms, "  After NRUN iter.:%14.5f   %13.5f   %13.5f\n", rmsrho, rmsz, rms);
	fprintf(filerms, "                  --------------------------------------------------------\n"); 
      }
      
      if(Nrunout != NULL) {
	sprintf(filename, "%s.txt", Nrunout);
	file = fopen(filename, "w");
	outpsi2xy(psi, file);
	fclose(file);

	sprintf(filename, "%s1d_rho.txt", Nrunout);
	file = fopen(filename, "w");
	outdenrho(psi, tmpz, file);
	fclose(file);     

	sprintf(filename, "%s1d_z.txt", Nrunout);
	file = fopen(filename, "w");
	outdenz(psi, tmprho, file);
	fclose(file);       
      }
   }
   
   if(rmsout != NULL) fclose(filerms);   
   
   fprintf(out, "                  --------------------------------------------------------\n\n");
   
   free_double_vector(rho);
   free_double_vector(z);

   free_double_vector(rho2);
   free_double_vector(z2);

   free_double_matrix(pot);
   free_double_matrix(psi);

   free_double_matrix(dpsirho);
   free_double_matrix(dpsiz);

   free_double_vector(calpharho);
   free_double_vector(calphaz);
   free_double_matrix(cbeta);
   free_double_vector(cgammarho);
   free_double_vector(cgammaz);

   free_double_matrix(tmprhoi);
   free_double_matrix(tmpzi);
   free_double_matrix(tmprhoj);
   free_double_matrix(tmpzj);
   
   free_double_vector(tmprho);
   free_double_vector(tmpz);
   
   clock_end = time(NULL);
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
   
   if((cfg_tmp = cfg_read("G0")) == NULL) {
     
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

      G0 = 4. * pi * as * Na * BOHR_RADIUS / aho;
   } else {
      G0 = atof(cfg_tmp);
   }

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
   
   if((cfg_tmp = cfg_read("NU")) == NULL) {
      fprintf(stderr, "NU is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   vnu = atof(cfg_tmp);
   
   if((cfg_tmp = cfg_read("LAMBDA")) == NULL) {
      fprintf(stderr, "LAMBDA is not defined in the configuration file.\n");
      exit(EXIT_FAILURE);
   }
   vlambda = atof(cfg_tmp);
   
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
   rmsout = cfg_read("RMSOUT");      
   initout = cfg_read("INITOUT");
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
   
   return;
}

/**
 *    Initialization of the space mesh, the potential, and the initial wave
 *    function.
 *    psi - array with the wave function values
 */
void init(double **psi) {
   long cnti, cntj;
   double vnu2, vlambda2;
   double pi3, cpsi;
   double tmp;
   
   if (opt == 1) par = 1.;
   else if (opt == 2) par = 2.;
   else{
      fprintf(stderr, "OPTION is not well defined in the configuration file\n");
      exit(EXIT_FAILURE);
   }
   
   vnu2 = vnu * vnu;
   vlambda2 = vlambda * vlambda;
   
   Nrho2 = Nrho / 2; Nz2 = Nz / 2;
   drho2 = drho * drho; dz2 = dz * dz;
   
   pi3 = pi * pi * pi;
   cpsi = sqrt(sqrt(pi3 / (vlambda * vnu2)));
   
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
            pot[cnti][cntj] = (vnu2 * rho2[cnti] + vlambda2 * z2[cntj]);
            tmp = exp(- 0.5 * (vnu * rho2[cnti] + vlambda * z2[cntj]));
            psi[cnti][cntj] = tmp / cpsi;
      }
   }
   
   return;
}

/**
 *    Crank-Nicolson scheme coefficients generation.
 */
void gencoef(void) {
   long cnti;
   
   Arho0 = 1. + dt / drho2;
   Az0 = 1. + dt / dz2;
   
   Arho0r = 1. - dt / drho2;
   Az0r = 1. - dt / dz2;
   
   Arho = - 0.5 * dt / drho2;
   Az = - 0.5 * dt / dz2;
   
   dArho = 0.25 * dt / drho;
   
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
void calcnorm(double *norm, double **psi, double **tmprho, double **tmpz) {
   int threadid;
   long cnti, cntj;
   
   #pragma omp parallel private(threadid, cnti, cntj)
   {
      threadid = omp_get_thread_num();
      
      #pragma omp for
      for(cnti = 0; cnti < Nrho; cnti ++) {
         for(cntj = 0; cntj < Nz; cntj ++) {
            tmpz[threadid][cntj] = psi[cnti][cntj] * psi[cnti][cntj];
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
void calcmuen(double *mu, double *en, double **psi, double **dpsirho, double **dpsiz, double **tmprhoi, double **tmpzi, double **tmprhoj, double **tmpzj) {
   int threadid;
   long cnti, cntj;
   double psi2, psi2lin, dpsi2;
   
   #pragma omp parallel private(threadid, cnti, cntj, psi2, psi2lin, dpsi2)
   {
      threadid = omp_get_thread_num();
      
      #pragma omp for
      for(cntj = 0; cntj < Nz; cntj ++) {
         for(cnti = 0; cnti < Nrho; cnti ++) {
            tmprhoi[threadid][cnti] = psi[cnti][cntj];
         }
         diff(drho, tmprhoi[threadid], tmprhoj[threadid], Nrho);
         for(cnti = 0; cnti < Nrho; cnti ++) {
            dpsirho[cnti][cntj] = tmprhoj[threadid][cnti];
         }
      }
      
      #pragma omp for
      for(cnti = 0; cnti < Nrho; cnti ++) {
         for(cntj = 0; cntj < Nz; cntj ++) {
            tmpzi[threadid][cntj] = psi[cnti][cntj];
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
            psi2 = psi[cnti][cntj] * psi[cnti][cntj];
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
void calcrms(double *rmsrho, double *rmsz, double *rms, double **psi, double **tmprhoi, double **tmpzi, double **tmprhoj, double **tmpzj) {
   int threadid;
   long cnti, cntj;
   double psi2;
   
   #pragma omp parallel private(threadid, cnti, cntj, psi2)
   {
      threadid = omp_get_thread_num();
      
      #pragma omp for
      for(cnti = 0; cnti < Nrho; cnti ++) {
         for(cntj = 0; cntj < Nz; cntj ++) {
            psi2 = psi[cnti][cntj] * psi[cnti][cntj];
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
void calcnu(double **psi) {
   long cnti, cntj;
   double psi2, psi2lin, tmp;
   
   #pragma omp parallel for private(cnti, cntj, psi2, psi2lin, tmp)
   for(cnti = 0; cnti < Nrho; cnti ++) {
      for(cntj = 0; cntj < Nz; cntj ++) {
         psi2 = psi[cnti][cntj] * psi[cnti][cntj];
         psi2lin = psi2 * G;
         tmp = dt * (pot[cnti][cntj] + psi2lin);
         psi[cnti][cntj] *= exp(- tmp);
      }
   }
   
   return;
}

/**
 *    Time propagation with respect to H2 (rho-part of the Laplacian).
 *    psi   - array with the wave function values
 *    cbeta - Crank-Nicolson scheme coefficients
 */
void calclurho(double **psi, double **cbeta) {
   int threadid;
   long cnti, cntj;
   double c;
   
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
void calcluz(double **psi, double **cbeta) {
   int threadid;
   long cnti, cntj;
   double c;
   
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

void outpsi2xy(double **psi, FILE *file) {
   long cnti, cntj;

   for(cnti = 0; cnti < Nrho; cnti += outstprho) {
      for(cntj = 0; cntj < Nz; cntj += outstpz) {
         fprintf(file, "%8le %8le %8le\n", rho[cnti], z[cntj], psi[cnti][cntj] * psi[cnti][cntj]);
      }
      fprintf(file, "\n");
      fflush(file);
   }
   return;
}

void outdenz(double **psi, double *tmp, FILE *file) {
   long cnti, cntj;

    for(cntj = 0; cntj < Nz; cntj += outstpz) {
        for(cnti = 0; cnti < Nrho; cnti ++) {
             tmp[cnti] = rho[cnti] * psi[cnti][cntj] * psi[cnti][cntj];
       }
       fprintf(file, "%8le %8le\n", z[cntj], 2. * pi * simpint(drho, tmp, Nrho));
       fflush(file);
    }
   return;    
}

void outdenrho(double **psi, double *tmp, FILE *file) {
   long cnti, cntj;  

    for(cnti = 0; cnti < Nrho; cnti += outstprho) {
       for(cntj = 0; cntj < Nz; cntj ++) {
             tmp[cntj] = psi[cnti][cntj] * psi[cnti][cntj];
       }
       fprintf(file, "%8le %8le\n", rho[cnti], simpint(dz, tmp, Nz));
       fflush(file);
    }
       return;
}
