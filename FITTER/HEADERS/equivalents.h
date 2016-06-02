#ifndef EQUIVALENTS_H
#define EQUIVALENTS_H

struct resampled **
ZN_average_slow( struct mom_info ***mominfoavg ,
		 struct input_params *INPARAMS ,
		 const struct mom_info **mominfo , 
		 const struct resampled **BOOTS ,
		 const int NSLICES ,
		 const bool actually_doit ,
		 const int Z_N ) ;

struct resampled *
momavg_single( struct mom_info **mominfoavg ,
	       int *NDATA ,
	       const struct mom_info *mominfo , 
	       const struct resampled *BOOTS ) ;

struct resampled **
momavg( struct mom_info ***mominfoavg ,
	struct input_params *INPARAMS ,
	const struct mom_info **mominfo , 
	const struct resampled **BOOTS ,
	const int NSLICES ,
	const bool actually_doit ) ;

// wrapper for momentum averaging
struct resampled **
momentum_average( struct mom_info ***mominfoavg ,
		  struct input_params *INPARAMS ,
		  const struct mom_info **mominfo , 
		  const struct resampled **BOOTS ,
		  const int NSLICES ,
		  const momavg_type type ) ;

#endif
