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
 *    Utility functions for integration (Simson's rule) and differentiation.
 */

#include "diffint.h"

/**
 *    Spatial 1D integration with Simpson's rule.
 *    h - space step
 *    f - array with the function values
 *    N - number of integration points
 */
double simpint(double h, double *f, long N) {
   long cnti;
   double sum, sumi, sumj, sumk;
   
   sumi = 0.; sumj = 0.; sumk = 0.;

   for(cnti = 1; cnti < N - 1; cnti += 2) {
      sumi += f[cnti];
      sumj += f[cnti - 1];
      sumk += f[cnti + 1];
   }
   
   sum = sumj + 4. * sumi + sumk;
   if(N % 2 == 0) sum += (5. * f[N - 1] + 8. * f[N - 2] - f[N - 3]) / 4.;

   return sum * h / 3.;
}

/**
 *    Richardson extrapolation formula for calculation of space
 *    derivatives.
 *    h  - space step
 *    f  - array with the function values
 *    df - array with the first derivatives of the function
 *    N  - number of space mesh points
 */
void diff(double h, double *f, double *df, long N) {
   long cnti;
   
   df[0] = 0.;
   df[1] = (f[2] - f[0]) / (2. * h);
   
   for(cnti = 2; cnti < N - 2; cnti ++) {
      df[cnti] = (f[cnti - 2] - 8. * f[cnti - 1] + 8. * f[cnti + 1] - f[cnti + 2]) / (12. * h); 
   }
   
   df[N - 2] = (f[N - 1] - f[N - 3]) / (2. * h);
   df[N - 1] = 0.;
   
   return;
}
