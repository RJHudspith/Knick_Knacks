#ifndef READ_TMOMENTS_H
#define READ_TMOMENTS_H

struct resampled**
read_tmoments( double ***X ,
	       struct mom_info ***mominfo ,
	       struct input_params *INPARAMS ,
	       int *NSLICES ) ;

#endif
