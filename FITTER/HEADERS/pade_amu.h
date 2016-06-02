/**
   @file pade_amu.h
   @brief compute amu from a pade fit
 */
#ifndef PADE_AMU_H
#define PADE_AMU_H

void
pade_amu_eval( double **xavg ,
	       struct resampled **bootavg ,
	       struct mom_info **mominfo ,
	       struct input_params *INPARAMS ,
	       const int NSLICES ,
	       const int LT ,
	       const bool renormalise ) ;

#endif
