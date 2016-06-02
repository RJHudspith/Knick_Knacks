#ifndef SPEED_OF_LIGHT_H
#define SPEED_OF_LIGHT_H

void
speed_of_light_eval( double **xavg ,
		   struct resampled **bootavg ,
		   struct mom_info **mominfo ,
		   struct input_params *INPARAMS ,
		   const int NSLICES ,
		   const int LT ) ;

#endif
