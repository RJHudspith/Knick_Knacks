#ifndef CORR_ELATION_H
#define CORR_ELATION_H

void
correlation( double **xavg ,
	     struct resampled **bootavg ,
	     struct mom_info **mominfo ,
	     struct input_params *INPARAMS ,
	     const int NSLICES ,
	     const int LT ) ;

#endif
