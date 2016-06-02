#ifndef TETRA_EVAL_H
#define TETRA_EVAL_H

void
tetra_eval( double **xavg ,
	    struct resampled **bootavg ,
	    struct mom_info **mominfo ,
	    struct mom_info *moms ,
	    struct input_params *INPARAMS ,
	    const int NSLICES ,
	    const int LT ) ;

#endif
