#ifndef FIT_AND_PLOT_H
#define FIT_AND_PLOT_H

struct resampled*
fit_data_plot_data( const struct resampled **BOOT ,
		    const double **X ,
		    const struct mom_info **mom_info ,
		    const struct input_params INPARAMS ,
		    const int NSLICES ,
		    const int LT ) ;

#endif
