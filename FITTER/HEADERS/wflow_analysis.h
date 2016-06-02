#ifndef WFLOW_ANALYSIS_H
#define WFLOW_ANALYSIS_H

void
analyse_w0( struct resampled **RAW , // need to resample this and root it ...
	    double *X ,
	    struct input_params INPARAMS ,// chiral data is in INPARAMS
	    const int NSLICES ) ;


#endif
