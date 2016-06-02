#ifndef J0CONTRIB
#define J0CONTRIB

void
j0_eval( double **xavg ,
	 struct resampled **bootavg ,
	 struct mom_info **mominfo ,
	 struct input_params *INPARAMS ,
	 const int NSLICES ,
	 const int LT ,
	 const bool renormalise ) ;

#endif
