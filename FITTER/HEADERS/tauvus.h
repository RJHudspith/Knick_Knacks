#ifndef TAUVUS_H
#define TAUVUS_H

void
tauvus_eval( double **xavg ,
	     struct resampled **bootavg ,
	     struct mom_info **mominfo ,
	     struct input_params *INPARAMS ,
	     const int NSLICES ,
	     const int LT ,
	     const bool renormalise ) ;

#endif
