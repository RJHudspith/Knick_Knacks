#ifndef READ_TETRAS_H
#define READ_TETRAS_H

struct resampled**
read_tetras( double ***X , // is the tiem
	     struct input_params *INPARAMS ,
	     struct mom_info **mom ,
	     int *NSLICES ,
	     const bool tfold ) ;

#endif
