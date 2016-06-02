#ifndef READ_GEVPCORR_H
#define READ_GEVPCORR_H

struct resampled**
read_gevpcorr( double ***X , // is the tiem
	       struct input_params *INPARAMS ,
	       struct mom_info **mom ,
	       int *NSLICES ,
	       const bool tfold ) ;

struct resampled**
read_gevpcorr2( double ***X , // is the tiem
		struct input_params *INPARAMS ,
		struct mom_info **mom ,
		int *NSLICES ,
		const bool tfold ) ;

#endif
