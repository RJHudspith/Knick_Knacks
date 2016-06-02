#ifndef READ_SIMCORRS_H
#define READ_SIMCORRS_H

struct resampled**
read_corr( double ***X , // is the momentum^2
	   struct input_params *INPARAMS ,
	   struct mom_info **mom ,
	   int *NSLICES ,
	   const int src ,
	   const int snk ,
	   const bool tfold ) ;

struct resampled**
read_simcorr( double ***X , // is the tiem
	      struct input_params *INPARAMS ,
	      struct mom_info **mom ,
	      int *NSLICES ,
	      const bool tfold ) ;

struct resampled**
read_simaxcorr( double ***X , // is the tiem
		struct input_params *INPARAMS ,
		struct mom_info **mom ,
		int *NSLICES ,
		const bool tfold ) ;

struct resampled**
read_simveccorr( double ***X , // is the tiem
		 struct input_params *INPARAMS ,
		 struct mom_info **mom ,
		 int *NSLICES ,
		 const bool tfold ) ;

struct resampled**
read_simbarcorr( double ***X , // is the tiem
		 struct input_params *INPARAMS ,
		 struct mom_info **mom ,
		 int *NSLICES ,
		 const bool tfold ) ;

struct resampled**
read_simvecbarcorr( double ***X , // is the tiem
		    struct input_params *INPARAMS ,
		    struct mom_info **mom ,
		    int *NSLICES ,
		    const bool tfold ) ;

#endif
