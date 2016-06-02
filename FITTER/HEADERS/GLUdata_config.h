#ifndef GLUDATA_CONFIG_H
#define GLUDATA_CONFIG_H

struct resampled**
read_GLUprop_config( struct mom_info ***mominfo ,
		     double ***X , // is the momentum^2
		     struct input_params *INPARAMS ,
		     int *NSLICES ,
		     const momtype mom_type ) ;

#endif
