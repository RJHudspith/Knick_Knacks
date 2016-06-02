#ifndef FLAVOUR_COMBINATION_H
#define FLAVOUR_COMBINATION_H

void
flavour_combination_eval( double **xavg ,
			  struct resampled **bootavg ,
			  struct mom_info **mominfo ,
			  struct input_params *INPARAMS ,
			  const int NSLICES ,
			  const int LT ,
			  const bool renormalise ) ;

#endif
