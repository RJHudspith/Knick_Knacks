#ifndef GLUFIT_CYLINDER_H
#define GLUFIT_CYLINDER_H

void
recylinder( struct input_params *INPARAMS ,    // change INPARAMS -> NDATA
	    struct resampled **boots ,
	    struct mom_info **mominfo ,    // pass by reference
	    double **x ,
	    const int NSLICES ,
	    const int LT ,
	    const int ND ,
	    const double *widths ) ;

#endif
