/**
   @file conformal_amu.h
   @brief conformal polynomial stuff
 */
#ifndef CONFORMAL_AMU_H
#define CONFORMAL_AMU_H

void
conformal_amu_eval( double **xavg ,
		    struct resampled **bootavg ,
		    struct mom_info **mominfo ,
		    struct input_params *INPARAMS ,
		    const int NSLICES ,
		    const int LT ,
		    const bool renormalise ) ;

#endif
