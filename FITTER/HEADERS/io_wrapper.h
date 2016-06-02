#ifndef IO_WRAPPER_H
#define IO_WRAPPER_H

int
read_data( double ***X ,
	   struct resampled ***RAW ,
	   struct mom_info ***mominfo ,
	   struct mom_info **moms ,
	   struct input_params *INPARAMS ,
	   int *NSLICES ,
	   const int LT ) ;

#endif
