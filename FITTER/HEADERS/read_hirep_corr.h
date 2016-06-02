#ifndef READ_HIREP_CORR
#define READ_HIREP_CORR

struct resampled**
read_hirep_corr( double ***X , // is the momentum^2
		 struct input_params *INPARAMS ,
		 struct mom_info **mom ,
		 int *NSLICES ,
		 const bool tfold ) ;

#endif
