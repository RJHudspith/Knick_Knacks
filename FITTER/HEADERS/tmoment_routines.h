#ifndef TMOMENT_ROUTINES_H
#define TMOMENT_ROUTINES_H

typedef enum {
  SYMMETRIC_MOMENTS ,
  PORTELLI_MOMENTS ,
  HPQCD_MOMENTS 
} moment_type ;

struct resampled*
compute_PIS( const struct resampled *tcorr ,
	     const int NMAX ,
	     const int LT , 
	     const moment_type moments ,
	     const int momtype ) ;

struct resampled**
compute_PADE( struct resampled **PIS ,
	      const int NSLICES ,
	      const int n ,
	      const int m ) ;

#endif
