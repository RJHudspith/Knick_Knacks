#ifndef GENERALISED_EVP_H
#define GENERALISED_EVP_H

void
gevp_eval( double **xavg ,
	   struct resampled **bootavg ,
	   struct mom_info **mominfo ,
	   struct input_params *INPARAMS ,
	   const int NSLICES ,
	   const int N , // correlation matrix length
	   const int LT ) ;

#endif
