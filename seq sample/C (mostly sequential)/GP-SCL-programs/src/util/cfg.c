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
 *    Utility functions for parsing of configuration files.
 */

#include "cfg.h"

/**
 *    Configuration file parsing.
 *    cfg_file - stands for a configuration file, which is supplied on
 *    a command line
 */
int cfg_init(char *cfg_file) {
   FILE *file;
   char buf[256];

   file = fopen(cfg_file, "r");
   if (! file) return 0;

   cfg_size = 0;
   while (fgets(buf, 256, file) != NULL) {
      if (sscanf(buf, "%s = %s", cfg_key[cfg_size], cfg_val[cfg_size]) == 2) cfg_size ++;
   }

   fclose(file);
   return cfg_size;
}

/**
 *    Configuration property value
 *    key - property
 */
char *cfg_read(char *key) {
   int i;

   for(i = 0; i < cfg_size; i ++)
      if (! strcmp(key, cfg_key[i])) return cfg_val[i];

   return NULL;
}
