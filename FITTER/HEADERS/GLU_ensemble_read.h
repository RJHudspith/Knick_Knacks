#ifndef GLU_ENSEMBLE_READ_H
#define GLU_ENSEMBLE_READ_H

struct resampled **
GLUdata( const char *filename , // filename (minus the number and .bin)
	 const int start ,      // trajectory start
	 const int end ,        // trajectory end
	 const int increment ,  // trajectory increment
	 const int NDATA   ,    // number of x's
	 const int NSLICES ,    // number of continguous arrays 
	 const bool POLYCORR ,   // are we measuring the POLY corr? 
	 const int ND           // dimensions of the MOM-list
	 ) ;

#endif
