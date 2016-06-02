#ifndef TIME_MOMENTS_H
#define TIME_MOMENTS_H

void
time_moment_eval( double **xavg ,
		  struct resampled **bootavg ,
		  struct mom_info **mominfo ,
		  struct input_params *INPARAMS ,
		  const int NSLICES ,
		  const int LT ,
		  const bool renormalise ) ;

#endif
