#ifndef DISPERSIONS_H
#define DISPERSIONS_H

void
dispersions_eval( double **xavg ,
		  struct resampled **bootavg ,
		  struct mom_info **mominfo ,
		  struct input_params *INPARAMS ,
		  const int NSLICES ,
		  const int LT ) ;

#endif
