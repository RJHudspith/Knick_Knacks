#ifndef B_CORRELATORS_H
#define B_CORRELATORS_H

void
bmeson_eval( double **xavg ,
	     struct resampled **bootavg ,
	     struct mom_info **mominfo ,
	     struct mom_info *moms ,
	     struct input_params *INPARAMS ,
	     const int NSLICES ,
	     const int LT ) ;

#endif
