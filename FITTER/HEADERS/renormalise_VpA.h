#ifndef RENORMALISE_VPA_H
#define RENORMALISE_VPA_H

void
renormalise_VpA( struct resampled **BOOT ,
		 const double **X , 
		 const struct inparams INPARAMS ) ;

void
renormalise_VmA( struct resampled **BOOT ,
		 const double **X , 
		 const struct inparams INPARAMS ) ;


#endif
