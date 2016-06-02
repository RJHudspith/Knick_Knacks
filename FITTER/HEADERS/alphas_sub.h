#ifndef ALPHAS_SUB
#define ALPHAS_SUB

void
alphas_sub( double **xavg ,
	    struct resampled **bootavg ,
	    struct mom_info **mominfo ,
	    struct input_params *INPARAMS ,
	    const int NSLICES ,
	    const int LT ,
	    const bool renormalise ) ;

#endif
