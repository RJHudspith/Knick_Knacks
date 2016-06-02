#ifndef MASS_EXTRAP_H
#define MASS_EXTRAP_H

void
mass_extrap( double **xavg ,
	     struct resampled **bootavg ,
	     struct mom_info **mominfo ,
	     struct input_params *INPARAMS ,
	     const int NSLICES ,
	     const int LT ) ;

#endif
