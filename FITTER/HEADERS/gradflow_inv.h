#ifndef GRADFLOW_INV_H
#define GRADFLOW_INV_H

int
gradflow_inv( double **X ,
	      struct resampled **BOOT ,
	      struct mom_info **mominfo ,
	      struct input_params *INPARAMS ,
	      const int NSLICES ,
	      const int LT ) ;

#endif
