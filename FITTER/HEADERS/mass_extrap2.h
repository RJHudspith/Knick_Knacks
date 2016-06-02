#ifndef MASS_EXTRAP2_H
#define MASS_EXTRAP2_H

void
mass_extrap2( double **xavg ,
	     struct resampled **bootavg ,
	     struct mom_info **mominfo ,
	     struct input_params *INPARAMS ,
	     const int NSLICES ,
	     const int LT ) ;

#endif
