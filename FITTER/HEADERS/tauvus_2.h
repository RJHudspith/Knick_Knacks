#ifndef TAUVUS2_H
#define TAUVUS2_H

void
tauvus2_eval( double **xavg ,
	      struct resampled **bootavg ,
	      struct mom_info **mominfo ,
	      struct input_params *INPARAMS ,
	      const int NSLICES ,
	      const int LT ,
	      const bool renormalise ) ;

#endif
