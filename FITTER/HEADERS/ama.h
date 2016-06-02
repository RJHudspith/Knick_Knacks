#ifndef AMA_H
#define AMA_H

void
ama_eval( double **xavg ,
	  struct resampled **bootavg ,
	  struct mom_info **mominfo ,
	  struct mom_info *moms ,
	  struct input_params *INPARAMS ,
	  const int NSLICES ,
	  const int LT ) ;

#endif
