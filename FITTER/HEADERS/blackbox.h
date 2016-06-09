#ifndef BLACKBOX_H
#define BLACKBOX_H

void
blackbox( const double *data ,
	  const int NDATA ,
	  const int NSTATES ,
	  double masses[ NSTATES ][ NDATA ] ) ;

struct resampled **
prony_effmass( const struct resampled **bootavg ,
	       const int *NDATA ,
	       const int NSLICES ,
	       const int NSTATES ,
	       const int STATE ) ;

#endif
