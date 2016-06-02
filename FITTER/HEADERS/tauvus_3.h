/**
   @file tauvus_3.h
   @brief time moment version of tau vus analysis
 */
#ifndef TAU_VUS_3_H
#define TAU_VUS_3_H

void
tauvus3_eval( double **xavg ,
	      struct resampled **bootavg ,
	      struct mom_info **mominfo ,
	      struct input_params *INPARAMS ,
	      const int NSLICES ,
	      const int LT ,
	      const bool renormalise ) ;

#endif
