#ifndef ANALYSIS_WRAPPER_H
#define ANALYSIS_WRAPPER_H

int
perform_analysis( double **xavg ,
		  struct resampled **BOOT , 
		  struct mom_info **mominfoavg ,
		  struct mom_info *moms ,
		  struct input_params *INPARAMS ,
		  const int NSLICES ,
		  const int LT ) ;

#endif
