/**
   @file statpot.c
   @brief static potential fit
 */
#ifndef STATPOT_H
#define STATPOT_H

void
statpot_eval( double **xavg ,
	      struct resampled **bootavg ,
	      struct mom_info **mominfo ,
	      const struct input_params *INPARAMS ,
	      const int NSLICES ,
	      const int LT ) ;

#endif
