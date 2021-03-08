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
#include <complex.h>

double *alloc_double_vector(long);
double complex *alloc_complex_vector(long);
double **alloc_double_matrix(long, long);
double ***alloc_double_tensor(long, long, long);
double complex ***alloc_complex_tensor(long, long, long);

void free_double_vector(double *);
void free_complex_vector(double complex *);
void free_double_matrix(double **);
void free_double_tensor(double ***);
void free_complex_tensor(double complex ***);
