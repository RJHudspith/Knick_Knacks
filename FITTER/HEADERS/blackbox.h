#ifndef BLACKBOX_H
#define BLACKBOX_H

struct resampled **
prony_effmass( const struct resampled **bootavg ,
	       const int *NDATA ,
	       const int NSLICES ,
	       const int NSTATES ,
	       const int STATE ) ;

#endif
