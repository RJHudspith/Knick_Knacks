#ifndef ARBITRARY_MOM_H
#define ARBITRARY_MOM_H

void
arbitrary_mom_eval( double **xavg ,
		    struct resampled **bootavg ,
		    struct mom_info **mominfo ,
		    struct input_params *INPARAMS ,
		    const int NSLICES ,
		    const int LT ,
		    const bool renormalise ) ;

#endif
