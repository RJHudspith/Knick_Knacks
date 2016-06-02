#ifndef ALPHAS_OPE_V1_H
#define ALPHAS_OPE_V1_H

void
alphas_v1_eval( double **xavg ,
		struct resampled **bootavg ,
		struct mom_info **mominfo ,
		struct input_params *INPARAMS ,
		const int NSLICES ,
		const int LT ,
		const bool renormalise ) ;

#endif
