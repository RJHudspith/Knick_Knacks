#ifndef BARYON_DEL_M
#define BARYON_DEL_M

void
baryon_mass_split_eval( double **xavg ,
			struct resampled **bootavg ,
			struct mom_info **mominfo ,
			struct mom_info *moms ,
			struct input_params *INPARAMS ,
			const int NSLICES ,
			const int LT ) ;

#endif
