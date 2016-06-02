#ifndef ALPHAS_OPE_H
#define ALPHAS_OPE_H

void
alphas_eval( double **xavg ,
	     struct resampled **bootavg ,
	     struct mom_info **mominfo ,
	     struct input_params *INPARAMS ,
	     const int NSLICES ,
	     const int LT ,
	     const bool renormalise ) ;

#endif
