#ifndef TETRA_EVAL_MESON_H
#define TETRA_EVAL_MESON_H

void
tetra_eval_meson( double **xavg ,
		  struct resampled **bootavg ,
		  struct mom_info **mominfo ,
		  struct mom_info *moms ,
		  struct input_params *INPARAMS ,
		  const int NSLICES ,
		  const int LT ) ;

void
tetra_eval_corr( double **xavg ,
		 struct resampled **bootavg ,
		 struct mom_info **mominfo ,
		 struct mom_info *moms ,
		 struct input_params *INPARAMS ,
		 const int NSLICES ,
		 const int LT ) ;

#endif
