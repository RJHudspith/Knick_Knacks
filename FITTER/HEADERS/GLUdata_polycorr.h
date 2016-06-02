#ifndef GLUDATA_POLYCORR_H
#define GLUDATA_POLYCORR_H

/**
   Reads the r^2 list for the polyakov loop correlator
 */
double*
xdata_double_polycorr( int *NDATA , 
		       int *NSLICES ,
		       const char *filename ,
		       const int start ) ;

#endif
