#ifndef RUN_AND_PLOT
#define RUN_AND_PLOT

void
run_fit_and_plot( struct resampled **RES ,
		  double **flatx ,
		  double **sigma ,
		  struct input_params INPARAMS ,
		  const int NSLICES ,
		  const int LT ) ;

#endif
