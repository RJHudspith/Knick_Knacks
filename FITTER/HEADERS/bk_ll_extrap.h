#ifndef BK_LL_EXTRAP_H
#define BK_LL_EXTRAP_H

void
ll_extrap_eval( double **xavg ,
		struct resampled **bootavg ,
		struct mom_info **mominfo ,
		struct mom_info *moms ,
		struct input_params *INPARAMS ,
		const int NSLICES ,
		const int LT ) ;

#endif
