#ifndef GENERALISED_EVP2_H
#define GENERALISED_EVP2_H

void
gevp2_eval( double **xavg ,
	    struct resampled **bootavg ,
	    struct mom_info **mominfo ,
	    struct input_params *INPARAMS ,
	    const int NSLICES ,
	    const int N , // correlation matrix length
	    const int LT ) ;

#endif
