#ifndef DFT_H
#define DFT_H

struct resampled **
momspace_corr( struct mom_info ***mom , 
	       double ***x ,
	       const struct resampled **tcorr ,
	       const struct input_params *INPARAMS ,
	       const int NSLICES , 
	       const int LT ) ;

#endif
