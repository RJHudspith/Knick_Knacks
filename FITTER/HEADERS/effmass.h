#ifndef EFFMASS_H
#define EFFMASS_H

struct resampled**
effective_mass( const struct resampled **bootavg ,
		const int *NDATA ,
		const int NSLICES ,
		const effmass_type type ) ;

#endif
